class FileReader:
    """Класс FileReader помогает читать из файла"""

    def __init__(self, path_to_file):
        self.path_to_file = path_to_file

    def read(self):
        try:
            with open(self.path_to_file, 'r') as f:
                return f.read()
        except IOError:
            return ''


if __name__ == '__main__':
    reader = FileReader('example.txt')
    print(reader.read())
