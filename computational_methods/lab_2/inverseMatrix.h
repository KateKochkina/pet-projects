#pragma once

#include "constants.h"
#include "qrMethod.h"


int computeInverseMatrix(const MyType* const* const A, MyType* const* const invA,
                         const int n)  // Вычисление обратой матрицы
{
    MyType** Q = new MyType*[n];  // Ортогональная матрица
    MyType** R = new MyType*[n];  // Верхнетреуголльная матрица
    for (int i = 0; i < n; ++i) {
        Q[i] = new MyType[n];
        R[i] = new MyType[n];
    }
    MyType* x_j = new MyType[n];   // j-й стобец искомой матрицы
    MyType* te_j = new MyType[n];  // j-й столбец единичной матрицы, умноженный слева на матрицу T

    qrDecomposition(A, Q, R, n);

    // T == tranposed(Q)
    // Решаем систему R * x_j = T * e_j
    for (int j = 0; j < n; ++j) {
        // Домножаем e_j слева на матрицу T
        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                MyType E_kj = (k == j) ? 1.0 : 0.0;
                sum += Q[k][i] * E_kj;
            }
            te_j[i] = sum;
        }

        backward((const MyType**) R, te_j, x_j, n);

        for (int i = 0; i < n; ++i) {
            invA[i][j] = x_j[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        delete[] R[i];
        delete[] Q[i];
    }
    delete[] Q;
    delete[] R;
    delete[] x_j;
    delete[] te_j;
}