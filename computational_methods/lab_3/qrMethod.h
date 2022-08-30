#pragma once

#include "constants.h"
#include "gaussMethod.h"
#include "inverseMatrix.h"
#include "matrix.h"


int qrMethod(const MyType* const* const A, const MyType* const b,
             MyType* const* const Q, MyType* const* const R,
             MyType* const x, const int n)  // Метод QR-разложения
{
    qrDecomposition(A, Q, R, n);

    // Если на диагонали есть нулевой элемент
    for (int k = 0; k < n; ++k) {
        if (fabs(R[k][k]) < EPS0) {
            return INVERTIBLE_MATRIX;
        }
    }

    // T == tranposed(Q)
    // Домножаем вектор правой части b слева на матрицу T
    MyType* tb = new MyType[n];

    for (int i = 0; i < n; ++i) {
        MyType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += Q[j][i] * b[j];
        }
        tb[i] = sum;
    }

    backward(R, tb, x, n);

    delete[] tb;

    return 0;
}


int qrMethod(const MyType* const* const A, const MyType* const b,
             MyType* const x, const int n)  // Метод QR-разложения
{
    MyType** Q = new MyType*[n];
    MyType** R = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        Q[i] = new MyType[n];
        R[i] = new MyType[n];
    }

    int res = qrMethod(A, b, Q, R, x, n);

    for (int i = 0; i < n; ++i) {
        delete[] Q[i];
        delete[] R[i];
    }
    delete[] Q;
    delete[] R;

    return res;
}