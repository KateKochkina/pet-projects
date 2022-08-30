#pragma once

#include "constants.h"
#include "gaussMethod.h"
#include "matrix.h"
#include "norm.h"
#include "vector.h"

const int VARIANT = 20;


void hessenberg(const MyType* const* const A, MyType* const* const H,
                const int n)  // Приведение матрицы к форме Хессенберга
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H[i][j] = A[i][j];
        }
    }

    for (int k = 1; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            if (fabs(H[i][k - 1]) < EPS0) {
                continue;  // Вращение не требуется
            }

            MyType denominator = sqrt(H[k][k - 1] * H[k][k - 1] + H[i][k - 1] * H[i][k - 1]);

            MyType a = H[k][k - 1] / denominator;
            MyType b = H[i][k - 1] / denominator;

            H[k][k - 1] = denominator;
            H[i][k - 1] = 0.0;
            for (int j = k; j < n; ++j) {
                MyType H_kj = H[k][j];
                MyType H_ij = H[i][j];
                H[k][j] = a * H_kj + b * H_ij;
                H[i][j] = -b * H_kj + a * H_ij;
            }

            H[k - 1][k] = denominator;
            H[k - 1][i] = 0.0;
            for (int j = k; j < n; ++j) {
                MyType H_jk = H[j][k];
                MyType H_ji = H[j][i];
                H[j][k] = a * H_jk + b * H_ji;
                H[j][i] = -b * H_jk + a * H_ji;
            }
        }
    }
}


int eigenValues(const MyType* const* const A, MyType* const eigVals, const int n,
                bool withShift = true, bool isHessenberg = true)  // Степенной метод (QR-разложения)
{
    MyType** A_copy = new MyType*[n];
    MyType** Q = new MyType*[n];
    MyType** R = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        Q[i] = new MyType[n];
        R[i] = new MyType[n];
        for (int j = 0; j < n; j++) {
            A_copy[i][j] = A[i][j];
        }
    }
    int* iter_count = new int[n];
    int* mult_count = new int[n];

    MyType max;
    for (int k = n - 1; k > 0; --k) {
        iter_count[k] = 0;
        mult_count[k] = 0;
        MyType eigVal_k = 0.0;

        do {
            ++iter_count[k];

            MyType d = A_copy[k][k];
            if (withShift) {
                for (int i = 0; i <= k; ++i) {
                    A_copy[i][i] -= d;
                    eigVal_k = eigVals[k];
                }
            }

            if (isHessenberg) {
                mult_count[k] += qrDecompositionHessenberg(A_copy, Q, R, k + 1);
            } else {
                mult_count[k] += qrDecomposition(A_copy, Q, R, k + 1);
            }

            // Если на диагонали есть нулевой элемент
            for (int i = 0; i < n; ++i) {
                if (fabs(R[i][i]) < EPS0) {
                    for (int j = 0; j < n; ++j) {
                        delete[] A_copy[j];
                        delete[] Q[j];
                        delete[] R[j];
                    }
                    delete[] A_copy;
                    delete[] Q;
                    delete[] R;
                    delete[] iter_count;
                    delete[] mult_count;

                    return INVERTIBLE_MATRIX;
                }
            }

            multiplyMatrix(R, Q, A_copy, k + 1);
            mult_count[k] += (k + 1) * (k + 1) * (k + 1);

            if (withShift) {
                for (int i = 0; i <= k; ++i) {
                    A_copy[i][i] += d;
                }
            }

            max = -INFINITY;
            for (int j = 0; j < k; ++j) {
                if (fabs(A_copy[k][j]) > max) {
                    max = fabs(A_copy[k][j]);
                }
            }
            if (iter_count[k] == 99 || iter_count[k] == 100) {
                cout << "k = " << k << ", итерация " << iter_count[k] << '\n';
                printMatrix(A_copy, k + 1);
            }
        } while (iter_count[k] < 100);//((max > sqrt(EPS0)) || (fabs(eigVal_k - eigVals[k]) > EPS));

        eigVals[k] = A_copy[k][k];
    }

    eigVals[0] = A_copy[0][0];

    cout << "Число итераций: ";
    cout << "{ " << iter_count[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << iter_count[i];
    }
    cout << " }\n";

    cout << "Число мултипликативных операций: ";
    cout << "{ " << mult_count[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << mult_count[i];
    }
    cout << " }\n";

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
        delete[] Q[i];
        delete[] R[i];
    }
    delete[] A_copy;
    delete[] Q;
    delete[] R;
    delete[] iter_count;
    delete[] mult_count;

    return 0;
}


int eigenVectors(const MyType* const* const A, const MyType* const eigVals,
                 MyType* const* const eigVecs, int n)  // Методом обратных итераций
{
    MyType** A_copy = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
    }
    MyType* x = new MyType[n];
    MyType* x_k = new MyType[n];
    MyType* y = new MyType[n];
    int* iter_count = new int[n];

    for (int i = 0; i < n; ++i) {
        iter_count[i] = 0;
        for (int j = 0; j < n; ++j) {
            x_k[j] = i == j ? 1.0 : 0.0;
        }

        for (int k = 0; k < n; ++k) {
            A_copy[k][k] = A[k][k] - eigVals[i];
        }

        do {
            if (iter_count[i] > 0) {
                for (int j = 0; j < n; ++j) {
                    x_k[j] = x[j];
                }
            }
            ++iter_count[i];

            int result = gaussMethod(A_copy, x_k, y, n);
            if (result != 0) {
                return result;
            }

            MyType y_norm = normVector(y, n);
            for (int j = 0; j < n; ++j) {
                x[j] = y[j] / y_norm;
            }
        } while ((1 - fabs(vectorDotProduct(x, x_k, n) / normVector(x, n) / normVector(x_k, n))) > EPS);

        for (int j = 0; j < n; ++j) {
            eigVecs[j][i] = x[j];
        }
    }

    cout << "Число итераций: ";
    cout << "{ " << iter_count[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << iter_count[i];
    }
    cout << " }\n";

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
    }
    delete[] A_copy;
    delete[] x;
    delete[] x_k;
    delete[] y;
    delete[] iter_count;

    return 0;
}


MyType rayleighRatio(const MyType* const* const A, MyType* const eigVals,
                     MyType* const* const eigVecs,
                     const int n)  // Комбинированный метод (соотношение Рэлея)
{
    MyType** A_copy = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
    }
    MyType* x = new MyType[n];
    MyType* x_k = new MyType[n];
    MyType* y = new MyType[n];
    MyType* temp = new MyType[n];
    int* iter_count = new int[n];

    for (int i = 0; i < n; ++i) {
        iter_count[i] = 0;
        MyType eigVal_i = 0.0;

        // Выбираем первые приближения собственных векторов
        if (VARIANT == 0) {  // 0 - тестовый пример
            switch (i) {
                case 0:
                    x_k[0] = 0.8;
                    x_k[1] = 0;
                    x_k[2] = 0.2;
                    x_k[3] = 0.4;
                    break;
                case 1:
                    x_k[0] = 0.0;
                    x_k[1] = -0.7;
                    x_k[2] = 0.6;
                    x_k[3] = -0.3;
                    break;
                case 2:
                    x_k[0] = -0.5;
                    x_k[1] = 0;
                    x_k[2] = 0.4;
                    x_k[3] = 0.7;
                    break;
                case 3:
                    x_k[0] = 0.0;
                    x_k[1] = -0.7;
                    x_k[2] = -0.6;
                    x_k[3] = 0.3;
                    break;
                default:
                    break;
            }
        } else if (VARIANT == 20) {
            switch (i) {
                case 0:
                    x_k[0] = 0.0;
                    x_k[1] = 0.0;
                    x_k[2] = 1.0;
                    x_k[3] = 0.0;
                    break;
                case 1:
                    x_k[0] = 0.9;
                    x_k[1] = 0.3;
                    x_k[2] = 0.0;
                    x_k[3] = 0.0;
                    break;
                case 2:
                    x_k[0] = 0.3;
                    x_k[1] = -0.9;
                    x_k[2] = 0.0;
                    x_k[3] = 0.2;
                    break;
                case 3:
                    x_k[0] = 0.0;
                    x_k[1] = 0.0;
                    x_k[2] = 0.0;
                    x_k[3] = 1.0;
                    break;
                default:
                    break;
            }
        }

        do {
            if (iter_count[i] > 0) {
                for (int j = 0; j < n; ++j) {
                    x_k[j] = x[j];
                }
                eigVal_i = eigVals[i];
            }
            ++iter_count[i];

            matrixMulByVector(A, x_k, temp, n);
            eigVals[i] = vectorDotProduct(temp, x_k, n);

            for (int k = 0; k < n; ++k) {
                A_copy[k][k] = A[k][k] - eigVals[i];
            }

            int result = gaussMethod(A_copy, x_k, y, n);
            if (result != 0) {
                return result;
            }

            MyType y_norm = normVector(y, n);
            for (int j = 0; j < n; ++j) {
                x[j] = y[j] / y_norm;
            }
        } while (((1 - fabs(vectorDotProduct(x, x_k, n) / normVector(x, n) / normVector(x_k, n))) > EPS) ||
                 (fabs(eigVal_i - eigVals[i]) > sqrt(EPS0)));

        for (int j = 0; j < n; ++j) {
            eigVecs[j][i] = x[j];
        }
    }

    cout << "Число итераций: ";
    cout << "{ " << iter_count[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << iter_count[i];
    }
    cout << " }\n";

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
    }
    delete[] A_copy;
    delete[] x;
    delete[] x_k;
    delete[] y;
    delete[] temp;
    delete[] iter_count;

    return 0;
}


void checkEig(const MyType* const* const A, const MyType* const eigVals,
              const MyType* const* const eigVecs, const int n)  // Проверка
{
    MyType** A_copy = new MyType*[n];
    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
    }

    cout << "\nПроверка:\n";
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            A_copy[k][k] = A[k][k] - eigVals[i];
        }
        cout << "det(A - lambda_" << i + 1 << " E) = "
             << fabs(computeMatrixDet(A_copy, n)) << '\n';
    }

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
    }
    delete[] A_copy;
}