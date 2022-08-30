#include "game.hpp"

#define PLAYER_DIED 1


void Player::put_on_equipment(const std::string &name) {
    equipment[name] = true;
    if (name == "armor") {
        ARM += 3;
        WGT += 3;
    } else if (name == "helmet") {
        ARM += 3;
        WGT += 2;
    } else if (name == "shield") {
        ARM += 5;
        WGT += 7;
    } else if (name == "pants") {
        ARM += 1;
        WGT += 1;
    } else if (name == "T-Shirt") {
        ARM += 1;
        WGT += 1;
    }
}

void Player::take_off_equipment(const std::string &name) {
    equipment[name] = false;
    if (name == "armor") {
        ARM -= 3;
        WGT -= 3;
    } else if (name == "helmet") {
        ARM -= 3;
        WGT -= 2;
    } else if (name == "shield") {
        ARM -= 5;
        WGT -= 7;
    } else if (name == "pants") {
        ARM -= 1;
        WGT -= 1;
    } else if (name == "T-Shirt") {
        ARM -= 1;
        WGT -= 1;
    }
}

void Events::create_map_from_file(std::ifstream &file) {
    file >> N >> M;
    create_map();

    std::string cell_content;
    size_t pos_X, pos_Y;

    while (file >> pos_X >> pos_Y >> cell_content) {
        game_map[pos_X][pos_Y] = cell_content;
    }
}

void Events::begin() {
    std::cout << "Supported actions:" << std::endl;

    if (N == 1 && M == 1) {
        std::cout << std::endl;
    } else {
        if (N != 1) {
            std::cout << " * move right" << std::endl;
        }
        if (M != 1) {
            std::cout << " * move up" << std::endl;
        }
    }

    std::cout << "0 x 0, hp: 100";
    if (stage_num == 2) {
        std::cout << ", armor: 0";
    }
    std::cout << " > ";
}

int Events::actions_process(Player *pl) {
    std::string act_1, act_2;

    while (std::cin >> act_1 >> act_2) {
        std::cout << std::endl;

        if (act_1 == "pick") {
            pl->put_on_equipment(act_2);
            std::cout << "clothes worn" << std::endl;

            game_map[pl->get_X()][pl->get_Y()] = "";
            actions_print(pl, "");
        }

        if (act_1 == "throw") {
            pl->take_off_equipment(act_2);
            std::cout << "the " << act_2 << " is thrown out" << std::endl;

            actions_print(pl, game_map[pl->get_X()][pl->get_Y()]);
        }

        if (act_1 == "move") {
            game_map[pl->get_X()][pl->get_Y()] = "";

            int dirX = 0, dirY = 0;

            if (act_2 == "left") {
                dirX = -1;
            } else if (act_2 == "right") {
                dirX = 1;
            } else if (act_2 == "down") {
                dirY = -1;
            } else if (act_2 == "up") {
                dirY = 1;
            }
            pl->move(dirX, dirY);

            std::string cell_content = game_map[pl->get_X()][pl->get_Y()];

            if (cell_content != "") {
                if (cell_content == "wolf") {
                    Wolf wolf;

                    if (battle_start(pl, &wolf, cell_content) == PLAYER_DIED) {
                        return PLAYER_DIED;
                    }
                } else if (cell_content == "dog") {
                    Dog dog;

                    if (battle_start(pl, &dog, cell_content) == PLAYER_DIED) {
                        return PLAYER_DIED;
                    }
                } else if (cell_content == "rat") {
                    Rat rat;

                    if (battle_start(pl, &rat, cell_content) == PLAYER_DIED) {
                        return PLAYER_DIED;
                    }
                } else {
                    std::cout << cell_content <<" found" << std::endl;
                    actions_print(pl, cell_content);
                }
            } else {
                std::cout << "moved" << std::endl;
                actions_print(pl, "");
            }
        }
    }

    return 0;
}

void Events::actions_print(Player *pl, const std::string &cell_content) {
    std::cout << "Supported actions:" << std::endl;

    size_t x = pl->get_X();
    size_t y = pl->get_Y();

    if (x != 0) {
        std::cout << " * move left" << std::endl;
    }
    if (x != N - 1) {
        std::cout << " * move right" << std::endl;
    }
    if (y != 0) {
        std::cout << " * move down" << std::endl;
    }
    if (y != M - 1) {
        std::cout << " * move up" << std::endl;
    }

    if (cell_content != "") {
        std::map<std::string, bool> equipment = pl->get_equipment();

        if (equipment[cell_content] == false) {
            int of_this_WGT = 1;

            if (cell_content == "armor") {
                of_this_WGT = 3;
            }
            if (cell_content == "helmet") {
                of_this_WGT = 2;
            }
            if (cell_content == "shield") {
                of_this_WGT = 7;
            }

            if (pl->get_WGT() + of_this_WGT <= pl->get_MAX_WGT()) {
                std::cout << " * pick " << cell_content << std::endl;
            }
        }

        for (auto it = equipment.begin(); it != equipment.end(); ++it) {
            if (it->second == true) {
                std::cout << " * throw " << it->first << std::endl;
            }
        }
    }

    std::cout << x << " x " << y << ", hp: " << pl->get_HP();
    if (stage_num == 2) {
        std::cout << ", armor: " << pl->get_ARM();
    }
    std::cout << " > ";
}

int Events::battle_start(Player *pl, Enemies *enemy, const std::string &cell_content) {
    std::cout << cell_content << " found, " << enemy->get_HP() << " hp" << std::endl;
    std::cout << "Supported actions:\n * kick enemy" << std::endl;

    std::cout << pl->get_X() << " x " << pl->get_Y() << ", hp: " << pl->get_HP();
    if (stage_num == 2) {
        std::cout << ", armor: " << pl->get_ARM();
    }
    std::cout << " > " << std::endl;

    std::string half_action;
    while (enemy->get_HP() > 0) {
        std::cin >> half_action >> half_action;

        if (battle_process(pl, enemy) == PLAYER_DIED) {
            return PLAYER_DIED;
        }
    }

    actions_print(pl, "");

    return 0;
}

int Events::battle_process(Player *pl, Enemies *enemy) {
    enemy->damaged_by(pl->get_DMG());

    if (enemy->get_HP() <= 0) {
        std::cout << "enemy killed" << std::endl;
    } else {
        if (enemy->get_DMG() > pl->get_ARM()) {
            pl->damaged_by(enemy->get_DMG() - pl->get_ARM());
        } else {
            pl->damaged_by(1);
        }

        if (pl->get_HP() <= 0) {
            std::cout << "player died" << std::endl;
            return PLAYER_DIED;
        }

        std::cout << "enemy kicked. Enemy hp: " << enemy->get_HP() << std::endl;
        std::cout << "Supported actions:\n * kick enemy" << std::endl;

        std::cout << pl->get_X() << " x " << pl->get_Y() << ", hp: " << pl->get_HP();
        if (stage_num == 2) {
            std::cout << ", armor: " << pl->get_ARM();
        }
        std::cout << " > " << std::endl;
    }

    return 0;
}
