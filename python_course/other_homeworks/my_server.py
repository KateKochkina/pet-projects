import asyncio

dict = {}

class ClientServerProtocol(asyncio.Protocol):
    def connection_made(self, transport):
        self.transport = transport

    def data_received(self, data):
        resp = process_data(data.decode())
        self.transport.write(resp.encode())


def process_data(data):
    list = data.split()

    if list[0] == 'put':
        key, value, timestamp = list[1], float(list[2]), int(list[3])
        if not key in dict:
            dict[key] = {}
        dict[key][timestamp] = value
        print(dict)
        return 'ok\n\n'

    if list[0] == 'get':
        recv = 'ok\n'
        key = list[1]
        if key == '*':
            for key in dict:
                for timestamp in dict[key]:
                    recv += '{} {} {}\n'.format(key, dict[key][timestamp], timestamp)
            recv += '\n'
            # print(recv)
            return recv
        if not key in dict:
            return 'ok\n\n'
        for timestamp in dict[key]:
            recv += '{} {} {}\n'.format(key, dict[key][timestamp], timestamp)
        recv += '\n'
        # print(recv)
        return recv

    return 'error\nwrong command\n\n'


def run_server(host, port):
    loop = asyncio.get_event_loop()
    coro = loop.create_server(
        ClientServerProtocol,
        host, port
    )

    server = loop.run_until_complete(coro)

    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass

    server.close()
    loop.run_until_complete(server.wait_closed())
    loop.close()


if __name__ == '__main__':
    run_server('127.0.0.1', 8888)
