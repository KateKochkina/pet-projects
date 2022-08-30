#pragma once

#include "constants.h"


MyType computeNorm1(const MyType* const x, const int n)  // Октаэдрическая норма вектора
{
    MyType sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += fabs(x[i]);
    }
    return sum;
}


MyType computeMatrixNorm1(const MyType* const* const A, const int n)  // Октаэдрическая норма матрицы
{
    MyType max = 0.0;
    for (int j = 0; j < n; ++j) {
        MyType sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += fabs(A[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}


MyType computeNorm2(const MyType* const x, const int n)  // Евклидова норма вектора
{
    MyType sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}


MyType computeNormInf(const MyType* const x, const int n)  // Кубическая норма вектора
{
    MyType max = 0.0;
    for (int i = 0; i < n; ++i) {
        if (fabs(x[i]) > max) {
            max = fabs(x[i]);
        }
    }
    return max;
}


MyType computeMatrixNormInf(const MyType* const* const A, const int n)  // Кубическая норма матрицы
{
    MyType max = 0.0;
    for (int i = 0; i < n; ++i) {
        MyType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += fabs(A[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}

#ifdef USE_NORM_1
MyType normVector(const MyType* const x, const int n) {
    return computeNorm1(x, n);
}
MyType normMatrix(const MyType* const* const A, const int n) {
    return computeMatrixNorm1(A, n);
}
#endif

#ifdef USE_NORM_2
MyType normVector(const MyType* const x, const int n) {
    return computeNorm2(x, n);
}
MyType normMatrix(const MyType* const* const A, const int n) {
    throw std::logic_error("Norm matrix 2 is not implemented.");
}
#endif

#ifdef USE_NORM_INF
MyType normVector(const MyType* const x, const int n) {
    return computeNormInf(x, n);
}
MyType normMatrix(const MyType* const* const A, const int n) {
    return computeMatrixNormInf(A, n);
}
#endif