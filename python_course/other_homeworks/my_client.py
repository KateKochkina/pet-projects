import socket
import time


class Client:
    def __init__(self, host, port, timeout = None):
        self.host = host
        self.port = port
        self.timeout = timeout

    def put(self, key, value, timestamp = None):
        with socket.create_connection((self.host, self.port), self.timeout) as sock:
            if timestamp is None:
                timestamp = str(int(time.time()))
            try:
                sock.send("put {} {} {}\n".format(key, float(value), timestamp).encode("utf8"))
            except socket.error:
                raise ClientError

            answer = sock.recv(1024).decode("utf8")
            if answer == "error\nwrong command\n\n":
                raise ClientError

    def get(self, key):
        dict = {}
        with socket.create_connection((self.host, self.port), self.timeout) as sock:
            try:
                sock.send("get {}\n".format(key).encode("utf8"))
            except socket.error:
                raise ClientError

            answer = sock.recv(1024).decode("utf8")
            if answer == "error\nwrong command\n\n":
                raise ClientError

            for line in answer.split("\n"):
                if line == "ok":
                    pass
                list = line.split()
                if len(list) == 3:
                    if not list[0] in dict:
                        dict[list[0]] = []
                    dict[list[0]].append((int(list[2]), float(list[1])))

        return dict

class ClientError(Exception):
    def __init__(self, text = ""):
        ClientError.txt = text


if __name__ == "__main__":
    client = Client("127.0.0.1", 8888, timeout=15)

    client.put("palm.cpu", 0.5, timestamp=1150864247)
    client.put("palm.cpu", 2.0, timestamp=1150864248)
    client.put("palm.cpu", 0.5, timestamp=1150864248)

    client.put("eardrum.cpu", 3, timestamp=1150864250)
    client.put("eardrum.cpu", 4, timestamp=1150864251)
    client.put("eardrum.memory", 4200000)

    print(client.get("*"))
