import os.path
import uuid

class File:
    def __init__(self, file_path):
        self.file_path = file_path
        self.current_position = 0

        if not os.path.exists(self.file_path):
            open(self.file_path, 'w').close()

    def read(self):
        with open(self.file_path, 'r') as f:
            return f.read()

    def write(self, content):
        with open(self.file_path, 'w') as f:
            f.write(content)

    def __add__(self, obj):
        new_path = os.path.join(os.path.dirname(self.file_path), str(uuid.uuid4().hex))
        new_file = type(self)(new_path)
        new_file.write(self.read() + obj.read())
        return new_file

    def __str__(self):
        return self.file_path

    def __iter__(self):
        return self

    def __next__(self):
        with open(self.file_path, 'r') as f:
            f.seek(self.current_position)

            line = f.readline()
            if not line:
                self.current_position = 0
                raise StopIteration('EOF')

            self.current_position = f.tell()
            return line


if __name__ == '__main__':
    a = File("tmp/python/1.txt")
    a.write("row1\n")
    a.write("row2\n")
    b = File("tmp/python/2.txt")
    b.write("row3\n")
    b.write("row4\n")
    c = a + b

    for i in a:
        print(i)
    print("----")
    for i in b:
        print(i)
    print("----")
    for i in c:
        print(i)
