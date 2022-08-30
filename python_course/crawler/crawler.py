import asyncio
from aiohttp import ClientSession
from aioelasticsearch import Elasticsearch
from collections import deque
from bs4 import BeautifulSoup
import re
import json
import time


loop = asyncio.get_event_loop()


class Crawler:
    def __init__(self, root_url, max_tasks, max_rps):
        self.root_url = root_url
        self.max_tasks = max_tasks
        self.max_rps = max_rps
        self.q = asyncio.Queue()
        self.q.put_nowait(root_url)
        self.seen_urls = set(root_url)
        self.es = Elasticsearch()
        self.timer = deque()

    async def crawl(self):
        self.session = ClientSession(loop=loop)
        workers = [asyncio.Task(self.work()) for _ in range(self.max_tasks)]
        await self.q.join()
        for w in workers:
            w.cancel()
        await self.session.close()
        await self.es.close()

    async def work(self):
        while True:
            url = await self.q.get()
            await self.fetch(url)
            self.q.task_done()

    async def fetch(self, url):
        await self.is_rps_exceeded()
        async with self.session.get(url) as response:
            html = await response.text(encoding='utf-8')
        await self.index_page(url, html)
        links = await self.parse_links(html)
        for link in links.difference(self.seen_urls):
            self.q.put_nowait(link)
            print(link)
        self.seen_urls.update(links)

    async def index_page(self, url, html):
        soup = BeautifulSoup(html, features='html.parser')
        [x.extract() for x in soup.find_all(['title', 'script', 'style',
                                             'meta'])]
        text = re.sub('<[^>]+>', '', str(soup))
        text = re.sub('(\s){2,}', ' ', text)
        content_url = json.dumps({
            'url': url,
            'content': text,
        })
        await self.es.index(index='crawling', doc_type='text', \
                            body=content_url)

    async def parse_links(self, html):
        links = set()
        soup = BeautifulSoup(html, features='html.parser')
        for link in soup.find_all('a'):
            href = link['href']
            if '#' in href:
                href = href.split('#', 1)[0]
            if '../' in href:
                href = href.split('../', 1)[1]
            if self.root_url in href:
                links.add(href)
            elif 'https://' in href or 'http://' in href:
                continue
            else:
                links.add(f'{self.root_url}/{href}')
        return links

    async def is_rps_exceeded(self):
        while True:
            now = time.perf_counter()
            while self.timer:
                if now - self.timer[0] > 1:
                    self.timer.popleft()
                else:
                    break
            if len(self.timer) < self.max_rps:
                break
            await asyncio.sleep(0.01)
        self.timer.append(time.perf_counter())


def main():
    crawler = Crawler('http://docs.python.org/', max_tasks=100, max_rps=10)
    loop.run_until_complete(crawler.crawl())


if __name__ == '__main__':
    main()
