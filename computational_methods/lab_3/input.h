#pragma once

#include "constants.h"

#include <fstream>
#include <iostream>

using std::ifstream;
using std::string;


int readSettings(const char* const fileName, size_t& n, string& matrixName,
                 string& outputName)  // Чтение параметров из файла
{
    ifstream file;
    file.open(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: file with settings is not open!\n";
        return 1;
    }

    file >> n;
    string str;
    getline(file, str);
    file >> matrixName;
    getline(file, str);
    file >> outputName;

    file.close();
    return 0;
}


int readData(MyType* const* const A, const string& matrixName)  // Чтение данных из файлов
{
    ifstream file;
    file.open(matrixName);
    if (!file.is_open()) {
        std::cerr << "Error: file with matrix is not open!\n";
        return 2;
    }

    int n;
    file >> n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }

    file.close();
    return 0;
}