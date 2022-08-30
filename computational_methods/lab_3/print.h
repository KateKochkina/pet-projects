#pragma once

#include "constants.h"
#include <iostream>
#include <iomanip>

using std::cout;


void printSystem(const MyType *const *const A, const MyType *const b, const int n)  // Вывод системы уравнений
{
    cout << "Система уравнений:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j > 0) {
                cout << ((A[i][j] < 0) ? " - " : " + ");
            } else {
                cout << ((A[i][j] < 0) ? "-" : " ");
            }
            cout << std::setiosflags(std::ios::left)
                 << std::setiosflags(std::ios::fixed)
                 << std::setprecision(3)
                 << std::setw(8) << fabs(A[i][j])
                 << std::setprecision(6)
                 << std::resetiosflags(std::ios::fixed)
                 << std::resetiosflags(std::ios::left)
                 << " * x" << j + 1 << "\t";
        }
        cout << " = " << b[i] << '\n';
    }
}


void printVector(const MyType *const x, const int n)  // Вывод вектора решения
{
    cout << "{ " << x[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << x[i];
    }
    cout << " }\n";
}


void printMatrix(const MyType *const *const A, const int n)  // Вывод матрицы
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << std::setw(15) << (fabs(A[i][j]) < EPS0 ? 0 : A[i][j]);
        }
        cout << '\n';
    }
}


void printResidual(const MyType *const *const A, const MyType *const b,
                   const MyType *const x, const int n,
                   MyType normVector(const MyType *const, const int))  // Вывод нормы невязки
{
    MyType *residual = new MyType[n];

    for (int i = 0; i < n; ++i) {
        MyType bi = 0;
        for (int j = 0; j < n; ++j) {
            bi += A[i][j] * x[j];
        }
        residual[i] = fabs(bi - b[i]);
    }

    cout << "Норма вектора невязки: " << normVector(residual, n) << '\n';

    delete[] residual;
}


void printError(const MyType *const exact, const MyType *const x, const int n,
                MyType normVector(const MyType *const, const int))  // Вывод нормы ошибки
{
    MyType *error = new MyType[n];

    for (int i = 0; i < n; ++i) {
        error[i] = fabs(exact[i] - x[i]);
    }

    cout << "Норма вектора ошибки: " << normVector(error, n) << '\n';

    delete[] error;
}