#pragma once

#include <typeinfo>


typedef double MyType;

const MyType EPS0 = (typeid(MyType).name()[0] == 'd') ? 1e-16 : 1e-6;
const MyType EPS = 1e-6;  // Задаваемая точность
const char *const SETTINGS_PATH = "settings.dat";


//#define USE_NORM_1
#define USE_NORM_2
//#define USE_NORM_INF


// Возвращаемые коды функций:

const int INVERTIBLE_MATRIX = -1;  // Матрица вырождена, система несовместна
