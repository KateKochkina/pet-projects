#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <cstring>
#include <map>

using std::string;
using std::map;

#define INVALID_ARGS_COUNT -1
#define ERR_FILE_OPEN -2


map <string, size_t> & get_data_from_file(std::ifstream &file) {
    map <string, size_t> *data = new map <string, size_t>;

    size_t id = 0;
    string name;
    size_t count = 0;

    while (file >> id >> name >> count) {
        std::transform(name.begin(), name.end(), name.begin(), (int (*)(int))std::tolower);
        (*data)[name] += count;
    }

    return *data;
}


void print_data(const map <string, size_t> &data) {
    size_t id = 0;

    for (auto it = data.begin(); it != data.end(); ++it, ++id) {
        std::cout << id << '|' << it->first << '|' << it->second << std::endl;
    }
}


int main(int argc, const char** argv) {
    if (argc != 2) {
        return INVALID_ARGS_COUNT;
    }

    const char *path_to_file = argv[1];

    std::ifstream file(path_to_file);
    if (!file.is_open()) {
        std::cout << "Eror file open" << std::endl;
        return ERR_FILE_OPEN;
    }

    map <string, size_t> data;

    data = get_data_from_file(file);

    print_data(data);

    data.clear();
    file.close();
    return 0;
}
