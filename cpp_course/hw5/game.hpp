#ifndef PROJECT_INCLUDE_GAME_HPP_
#define PROJECT_INCLUDE_GAME_HPP_

#include <iostream>
#include <fstream>
#include <map>


class Player {
 public:
    Player() : X(0), Y(0), HP(100), DMG(1), ARM(0), WGT(0), MAX_WGT(20) {
        equipment["armor"] = false;
        equipment["helmt"] = false;
        equipment["shield"] = false;
        equipment["pants"] = false;
        equipment["T-Shirt"] = false;
    }

    int get_X() {
        return X;
    }
    int get_Y() {
        return Y;
    }
    int get_HP() {
        return HP;
    }
    int get_DMG() {
        return DMG;
    }
    int get_ARM() {
        return ARM;
    }
    int get_WGT() {
        return WGT;
    }
    int get_MAX_WGT() {
        return MAX_WGT;
    }
    const std::map<std::string, bool> & get_equipment() {
        return equipment;
    }

    void move(int dirX, int dirY) {
        X += dirX;
        Y += dirY;
    }

    void damaged_by(int val) {
        HP -= val;
    }

    void put_on_equipment(const std::string &name);

    void take_off_equipment(const std::string &name);

 private:
    size_t X, Y;
    int HP;
    int DMG;
    int ARM;
    int WGT;
    const int MAX_WGT;
    std::map<std::string, bool> equipment;
};


class Enemies {
 public:
    int get_HP() {
        return HP;
    }

    int get_DMG() {
        return base_DMG + weapon_DMG;
    }

    void damaged_by(int val) {
        HP -= val;
    }

 protected:
    int HP;
    int base_DMG;
    int weapon_DMG;
};


class Wolf : public Enemies {
 public:
    Wolf() {
        HP = 6;
        base_DMG = 10;
        weapon_DMG = 1;
    }
};

class Dog : public Enemies {
 public:
    Dog() {
        HP = 3;
        base_DMG = 5;
        weapon_DMG = 1;
    }
};

class Rat : public Enemies {
 public:
    Rat() {
        HP = 2;
        base_DMG = 3;
        weapon_DMG = 1;
    }
};


class Events {
 public:
    explicit Events(int stage_num) : stage_num(stage_num), N(0), M(0), game_map(NULL) {}

    ~Events() {
        for (size_t i = 0; i < N; ++i) {
            delete []game_map[i];
        }
        delete []game_map;
    }

    void create_map_from_file(std::ifstream &file);

    void begin();

    void play() {
        begin();
        Player player;
        actions_process(&player);
    }

 private:
    int stage_num;
    size_t N, M;
    std::string** game_map;

    void create_map() {
        game_map = new std::string*[N];
        for (size_t i = 0; i < N; ++i) {
            game_map[i] = new std::string[M];
        }
    }

    int actions_process(Player *pl);

    void actions_print(Player *pl, const std::string &cell_content);

    int battle_start(Player *pl, Enemies *enemy, const std::string &cell_content);

    int battle_process(Player *pl, Enemies *enemy);
};


#endif  // PROJECT_INCLUDE_GAME_HPP_
