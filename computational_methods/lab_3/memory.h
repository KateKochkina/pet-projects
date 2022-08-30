#pragma once

#include "constants.h"


void memoryAllocation(MyType**& A, MyType**& H, MyType*& eigVals, MyType**& eigVecs,
                      const int n)  // Выделение памяти
{
    A = new MyType*[n];
    H = new MyType*[n];
    eigVals = new MyType[n];
    eigVecs = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new MyType[n];
        H[i] = new MyType[n];
        eigVecs[i] = new MyType[n];
    }
}


void memoryDeallocation(MyType**& A, MyType**& H, MyType*& eigVals, MyType**& eigVecs,
                        const int n)  // Освобождение памяти
{
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] H[i];
        delete[] eigVecs[i];
    }
    delete[] A;
    delete[] H;
    delete[] eigVals;
    delete[] eigVecs;
}
