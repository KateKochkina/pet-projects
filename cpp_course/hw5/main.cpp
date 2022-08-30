#include <iostream>
#include <cstring>

#include "game.hpp"

#define INVALID_ARGS_COUNT 1
#define INVALID_MAP_KEY 2
#define ERROR_FILE_OPEN 3


int main(int argc, const char** argv) {
    if (argc != 3 && argc != 4) {
        return INVALID_ARGS_COUNT;
    }

    int stage_num = 1;
    const char *path_to_map_file = "\0";

    if (!strcmp(argv[1], "--map")) {
        path_to_map_file = argv[2];

        if (argc == 4 && !strcmp(argv[3], "--view-armor")) {
            stage_num = 2;
        }
    } else if (!strcmp(argv[1], "--view-armor")) {
        stage_num = 2;

        if (!strcmp(argv[2], "--map")) {
            path_to_map_file = argv[3];
        } else {
            return INVALID_MAP_KEY;
        }
    }

    std::ifstream file(path_to_map_file);
    if (!file.is_open()) {
        return ERROR_FILE_OPEN;
    }

    Events event(stage_num);

    event.create_map_from_file(file);
    file.close();

    event.play();

    return 0;
}
