import asyncio
import aiohttp
import sys

async def main():
    async with aiohttp.ClientSession() as session:
        params = {
            'q': ' '.join(sys.argv[1:]),
            'limit': 100,
            'offset': 0,
        }
        async with session.get('http://localhost:8080/api/v1/search', \
                               params=params) as response:
            print(await response.text())

loop = asyncio.get_event_loop()
loop.run_until_complete(main())
