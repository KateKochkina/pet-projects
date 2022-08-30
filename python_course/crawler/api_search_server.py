from aiohttp import web
from aioelasticsearch import Elasticsearch


async def search(request):
    query = request.rel_url.query
    body = {
        "size": query['limit'],
        "from": query['offset'],
        "query": {
            "match_phrase": {
                "content": query['q']
            }
        }
    }
    res = await es.search(index='crawling', doc_type='text', body=body)
    urls = []
    for hit in res['hits']['hits']:
        url = hit['_source']['url']
        urls.append(url)
    return web.Response(text='\n'.join(urls))

es = Elasticsearch()
app = web.Application()
app.add_routes([web.get('/api/v1/search', search)])

web.run_app(app)
