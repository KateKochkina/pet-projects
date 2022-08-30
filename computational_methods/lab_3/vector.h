#pragma once

#include "constants.h"


MyType vectorDotProduct(const MyType* const a, const MyType* const b,
                        const int n)  // Евклидово скалярное произведение
{
    MyType sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}


void vectorScalarMul(const MyType* const a, const MyType x, MyType* const result,
                     const int n)  // Умножение вектора на число
{
    for (int i = 0; i < n; ++i) {
        result[i] = a[i] * x;
    }
}

MyType* vectorAbsDiff(const MyType* const a, const MyType* const b, MyType* const result,
                   const int n)  // Разность двух векторов
{
    MyType* t;
    t = new MyType[n];

    for (int i = 0; i < n; ++i) {
        result[i] = fabs(a[i] - b[i]);
    }

    return t;
}
