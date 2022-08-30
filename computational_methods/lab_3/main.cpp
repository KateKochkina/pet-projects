#include <cmath>
#include <iostream>

#include "constants.h"
#include "eigen.h"
#include "input.h"
#include "memory.h"
#include "print.h"


int main() {
    setlocale(LC_ALL, "Russian");

    size_t n;  // Размер матрицы
    string matrixFileName, vectorFileName, outputFileName;
    int result = readSettings(SETTINGS_PATH, n, matrixFileName, outputFileName);
    if (result != 0) {
        return result;
    }

    MyType **A,         // Матрица коэффициентов
            **H,        // Матрица Хессенберга
            *eigVals,   // Массив собственных значений
            **eigVecs;  // Матрица, столбцами которой являются собственные векторы

    memoryAllocation(A, H, eigVals, eigVecs, n);

    result = readData(A, matrixFileName);
    if (result != 0) {
        return result;
    }

    cout << "Исходная матрица:\n";
    printMatrix(A, n);

    hessenberg(A, H, n);
    cout << "Матрица Хессенберга:\n";
    printMatrix(A, n);

    cout << "\nПоиск собственных значений методом QR-разложения:\n";
    result = eigenValues(A, eigVals, n, false, false);
    if (result == 0) {
        cout << "Приближенные собственные значения: ";
        printVector(eigVals, n);
    } else if (result == INVERTIBLE_MATRIX) {
        cout << "Матрица A вырождена, система несовместна.\n";
    } else {
        cout << "Unknown result: " << result << '\n';
    }

    cout << "\nПоиск собственных векторов методом обратных итераций:\n";
    result = eigenVectors(A, eigVals, eigVecs, n);
    if (result == 0) {
        cout << "Приближенные собственные векторы:\n";
        printMatrix(eigVecs, n);
    } else if (result == INVERTIBLE_MATRIX) {
        cout << "Матрица A вырождена, система несовместна.\n";
    } else {
        cout << "Unknown result: " << result << '\n';
    }

    checkEig(A, eigVals, eigVecs, n);

    cout << "\nПоиск собственных значений и собственных векторов с помощью отношения Рэлея:\n";
    result = rayleighRatio(A, eigVals, eigVecs, n);
    if (result == 0) {
        cout << "Приближенные собственные значения: ";
        printVector(eigVals, n);
        cout << "Приближенные собственные векторы:\n";
        printMatrix(eigVecs, n);
    } else if (result == INVERTIBLE_MATRIX) {
        cout << "Матрица A вырождена, система несовместна.\n";
    } else {
        cout << "Unknown result: " << result << '\n';
    }

    checkEig(A, eigVals, eigVecs, n);

    memoryDeallocation(A, H, eigVals, eigVecs, n);
    return 0;
}