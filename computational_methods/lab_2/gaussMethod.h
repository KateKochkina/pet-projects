#pragma once

#include "constants.h"


void backward(const MyType* const* const A, const MyType* const b,
              MyType* const x, const int n)  // Обратный ход метода Гаусса
{
    for (int i = 1; i < n; ++i) {
        x[i] = 0.0;
    }

    MyType t;
    for (int i = n - 1; i >= 0; --i) {
        t = b[i];
        for (int j = i + 1; j < n; ++j) {
            t -= A[i][j] * x[j];
        }

        x[i] = t / A[i][i];
        if (fabs(x[i]) < EPS0) {
            x[i] = 0.0;
        }
    }
}


int gaussMethod(const MyType* const* const A, const MyType* const b, MyType* const x,
                const int n)  // Метода Гаусса
{
    MyType** A_copy = new MyType*[n];
    MyType* b_copy = new MyType[n];

    // Копируем данные системы
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
        b_copy[i] = b[i];
    }

    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n; ++k) {
        // Частичный поиск ведущего элемента по столбцу
        int iMax = k;
        for (int i = k; i < n; ++i) {
            if (fabs(A_copy[i][k]) > fabs(A_copy[iMax][k])) {
                iMax = k;
            }
        }

        // Если на диагонали нулевой элемент
        if (fabs(A_copy[iMax][k]) < EPS0) {
            for (int i = 0; i < n; ++i) {
                delete[] A_copy[i];
            }
            delete[] A_copy;
            delete[] b_copy;

            return INVERTIBLE_MATRIX;
        }

        // Меняем строки
        if (iMax != k) {
            MyType* temp = A_copy[k];
            A_copy[k] = A_copy[iMax];
            A_copy[iMax] = temp;

            MyType t = b_copy[k];
            b_copy[k] = b_copy[iMax];
            b_copy[iMax] = t;
        }

        // Прямой ход
        for (int i = k + 1; i < n; ++i) {
            MyType t = A_copy[i][k] / A_copy[k][k];
            for (int j = k; j < n; ++j) {
                A_copy[i][j] -= A_copy[k][j] * t;
            }
            b_copy[i] -= b_copy[k] * t;
        }
    }

    backward(A_copy, b_copy, x, n);

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
    }
    delete[] A_copy;
    delete[] b_copy;

    return 0;
}