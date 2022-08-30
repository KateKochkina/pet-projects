#pragma once

#include <typeinfo>


typedef double MyType;

const MyType EPS0 = (typeid(MyType).name()[0] == 'd') ? 1e-12 : 1e-6;
const MyType EPS = 1e-8;  // Задаваемая точность

#define USE_NORM_INF  // иначе норма 1

// Возвращаемые коды функций
const int INVERTIBLE_MATRIX = -1;  // Матрица вырождена, система несовместна
