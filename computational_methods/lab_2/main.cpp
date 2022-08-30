#include "constants.h"
#include "inverseMatrix.h"
#include "norm.h"
#include "print.h"

#include <cmath>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

const char *const settingsFileName = "settings.dat";

int readSettings(const char *fileName, size_t &n, string &matrixName, string &vectorName,
                 string &outputName);

void memoryAllocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                      MyType *&x, int n);

void memoryDeallocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                        MyType *&x, int n);

int readData(MyType *const *A, MyType *b, const string &matrixName, const string &vectorName);

void correction(MyType *const *A, MyType *b, int n);

void LDU(const MyType *const *A, MyType *const *L, MyType *const *D, MyType *const *U, int n);

int fixedPointIterMethod(const MyType *const *A, const MyType *b, MyType *x, MyType tau, int n,
                         const MyType *exact);

int jacobiMethod(const MyType *const *A, const MyType *b, MyType *x, int n, const MyType *exact);

int seidelMethod(const MyType *const *A, MyType *const *L, MyType *const *D, MyType *const *U,
                 MyType *b, MyType *x, MyType omega, int n, const MyType *exact);

int seidelThreeDiag(const MyType *a, const MyType *b, const MyType *c, const MyType *d,
                    MyType *x, MyType omega, int n);


int main()  // Точка входа в приложение
{
    setlocale(LC_ALL, "Russian");

    size_t n;  // Размер матрицы
    string matrixFileName, vectorFileName, outputFileName;
    int result = readSettings(settingsFileName, n, matrixFileName,
                              vectorFileName, outputFileName);
    if (result != 0) {
        return result;
    }

    MyType **A,   // Матрица коэффициентов
            **L,  // Нижний треугольник матрицы А
            **D,  // Диагональ матрицы А
            **U,  // Верхний треугольник матрицы A
            *q,   // Вектор правой части
            *x;   // Вектор неизвестных

    memoryAllocation(A, L, D, U, q, x, n);

    result = readData(A, q, matrixFileName, vectorFileName);
    if (result != 0) {
        return result;
    }

    cout << "\nВариант 19, система " << matrixFileName[11] << ", тип "
         << ((typeid(MyType).name()[0] == 'd') ? "double" : "float")
         << ", кубическая норма, точность " << EPS << "\n\n";

    correction(A, q, n);
    printSystem(A, q, n);

    // Точное решение системы 2 варианта 19
    const MyType exact[4] = {6., -9., 12., 8.};

    // Решение системы методом простой итерации
    cout << "\nМетод простой итерации:" << '\n';
    MyType tau = 0.01384;
    cout << "tau = " << tau << " (оптимальный для данной системы)\n";
    fixedPointIterMethod(A, q, x, tau, n, exact);

    // Решение системы методом Якоби
    cout << "\nМетод Якоби:\n";
    jacobiMethod(A, q, x, n, exact);

    // Решение системы методом Зейделя
    cout << "\nМетод Зейделя:\n";
    LDU(A, L, D, U, n);
    seidelMethod(A, L, D, U, q, x, 1, n, exact);

    // Решение системы методом релаксации
    cout << "\nМетод релаксации:\n";
    MyType omega = 1.0;
    cout << "omega = " << omega << " (оптимальный для данной системы)\n";
    //seidelMethod(A, L, D, U, q, x, omega, n, exact);
    cout << "Совпадает с решением, полученным методом Зейделя\n";

    memoryDeallocation(A, L, D, U, q, x, n);

    // Решение методом релаксации для трехдиагональных матриц
    const int N_VARIANT = 19;
    const int m = 200 + N_VARIANT;

    MyType *a,   // Поддиагональ большой матрицы
            *b,  // Диагональ большой матрицы
            *c,  // Наддиагональ большой матрицы
            *d;  // Правая часть системы

    a = new MyType[m];
    b = new MyType[m];
    c = new MyType[m];
    d = new MyType[m];
    x = new MyType[m];

    cout << "\nМетод релаксации для трехдиагональной матрицы размера " << m << "x" << m << ":\n"
         << "a[i] = c[i] = 1; b[i] = 4; i = 1, 2, ..., m;\n"
            "d[1] = 6; d[i] = 10−2·(i mod 2), i = 2, 3, ..., m−1; d[m] = 9−3(n mod 2).\n";
    for (int i = 0; i < m; ++i) {
        a[i] = c[i] = 1;
        b[i] = 4;
        d[i] = 10 - 2 * ((i + 1) % 2);
    }
    d[0] = 6;
    d[m - 1] = 9 - 3 * (m % 2);
    seidelThreeDiag(a, b, c, d, x, 1, m);

    cout << "\nМетод релаксации для трехдиагональной матрицы размера " << m << "x" << m << ":\n"
         << "a[i] = 6,\tb[i] = 8,\tc[i] = 1,\td[i] = i\n";
    for (int i = 0; i < m; ++i) {
        a[i] = 6;
        b[i] = 8;
        c[i] = 1;
        d[i] = i;
    }
    seidelThreeDiag(a, b, c, d, x, 1, m);

    cout << "\nМетод релаксации для трехдиагональной матрицы размера " << m << "x" << m << ":\n"
         << "a[i] = 1,\tb[i] = 8,\tc[i] = 6,\td[i] = i\n";
    for (int i = 0; i < m; ++i) {
        a[i] = 1;
        b[i] = 8;
        c[i] = 6;
        d[i] = i;
    }
    seidelThreeDiag(a, b, c, d, x, 1, m);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] x;

    return 0;
}


int readSettings(const char *const fileName, size_t &n, string &matrixName,
                 string &vectorName, string &outputName)  // Чтение параметров из файла
{
    ifstream file;
    file.open(fileName);
    if (!file.is_open()) {
        cerr << "Error: file with settings is not open!\n";
        return 1;
    }

    file >> n;
    string str;
    getline(file, str);
    file >> matrixName;
    getline(file, str);
    file >> vectorName;
    getline(file, str);
    file >> outputName;

    file.close();
    return 0;
}


void memoryAllocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                      MyType *&x, const int n)  // Выделение памяти
{
    A = new MyType *[n];
    L = new MyType *[n];
    D = new MyType *[n];
    U = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new MyType[n];
        L[i] = new MyType[n];
        D[i] = new MyType[n];
        U[i] = new MyType[n];
    }
    b = new MyType[n];
    x = new MyType[n];
}


void memoryDeallocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                        MyType *&x, const int n)  // Освобождение памяти
{
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] L[i];
        delete[] D[i];
        delete[] U[i];
    }
    delete[] A;
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] b;
    delete[] x;
}


int readData(MyType *const *const A, MyType *const b, const string &matrixName,
             const string &vectorName)  // Чтение данных из файлов
{
    ifstream file;
    file.open(matrixName);
    if (!file.is_open()) {
        cerr << "Error: file with matrix is not open!\n";
        return 2;
    }

    int n;
    file >> n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }

    file.close();
    file.open(vectorName);
    if (!file.is_open()) {
        cerr << "Error: file with vector is not open!\n";
        return 3;
    }

    file >> n;
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }

    file.close();
    return 0;
}


void correction(MyType *const *const A, MyType *const b, const int n)  // Корректировка системы
{
    MyType *sum = new MyType[n];
    for (int i = 0; i < n; ++i) {
        sum[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                sum[i] += fabs(A[i][j]);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (A[i][i] < 0.0 && fabs(A[i][i]) > sum[i]) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = -A[i][j];
            }
            b[i] = -b[i];
        }
    }
    delete[] sum;
}


void LDU(const MyType *const *const A, MyType *const *const L, MyType *const *const D,
         MyType *const *const U, const int n)  // LDU-разложение матрицы
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                D[i][j] = A[i][j];
                L[i][j] = 0;
                U[i][j] = 0;
            } else if (i > j) {
                D[i][j] = 0;
                L[i][j] = A[i][j];
                U[i][j] = 0;
            } else {
                D[i][j] = 0;
                L[i][j] = 0;
                U[i][j] = A[i][j];
            }
        }
    }
}


int fixedPointIterMethod(const MyType *const *const A, const MyType *const b, MyType *const x,
                         const MyType tau, const int n,
                         const MyType *const exact)  // Метод простой итерации
{
    MyType **C,
            *y,
            *x_k,
            *diff;

    C = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        C[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            MyType e_ij = (i == j) ? 1 : 0;
            C[i][j] = e_ij - tau * A[i][j];
        }
        y[i] = tau * b[i];
        x_k[i] = y[i];
    }

    cout << "Вектор y: ";
    printVector(y, n);
    cout << "Матрица C:" << '\n';
    printMatrix(C, n);
    MyType cNorm = normMatrix(C, n);
    cout << "Норма матрицы C: " << cNorm << '\n';

    double iter_est = 0.0;
    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += C[i][j] * x_k[j];
            }
            x[i] = t + y[i];
        }

        //        {
        //            for (int i = 0; i < n; ++i) {
        //                diff[i] = x[i] - exact[i];
        //            }
        //
        //            ofstream file("test/output2.dat", std::ios::app);
        //            if (!file.is_open()) {
        //                cerr << "Error: file with matrix is not open!\n";
        //                return 1;
        //            }
        //            file << iter << ' ' << normVector(diff, n) << '\n';
        //        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
        if (iter == 1) {
            iter_est = log(EPS / normVector(diff, n) * (1 - cNorm)) / log(cNorm);
        }
    } while (normVector(diff, n) > (EPS * (1 - cNorm)) / cNorm);

    cout << "Оценка числа итераций: > " << iter_est << '\n';
    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    printError(exact, x, n, normVector);
    printResidual(A, b, x, n, normVector);

    for (int i = 0; i < n; ++i) {
        delete[] C[i];
    }
    delete[] C;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int jacobiMethod(const MyType *const *const A, const MyType *const b, MyType *const x,
                 const int n, const MyType *const exact)  // Метод Якоби
{
    MyType **C,
            *y,
            *x_k,
            *diff;

    C = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        C[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = (i == j) ? 0 : -A[i][j] / A[i][i];
        }
        y[i] = b[i] / A[i][i];
        x_k[i] = y[i];
    }

    cout << "Вектор y: ";
    printVector(y, n);
    cout << "Матрица C: " << '\n';
    printMatrix(C, n);
    MyType cNorm = normMatrix(C, n);
    cout << "Норма матрицы C: " << cNorm << '\n';

    double iter_est = 0.0;
    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += C[i][j] * x_k[j];
            }
            x[i] = t + y[i];
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
        if (iter == 1) {
            iter_est = log(EPS / normVector(diff, n) * (1 - cNorm)) / log(cNorm);
        }
    } while (normVector(diff, n) > (EPS * (1 - cNorm)) / cNorm);

    cout << "Оценка числа итераций: > " << iter_est << '\n';
    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    printError(exact, x, n, normVector);
    printResidual(A, b, x, n, normVector);

    for (int i = 0; i < n; ++i) {
        delete[] C[i];
    }
    delete[] C;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int seidelMethod(const MyType *const *const A, MyType *const *const L, MyType *const *const D,
                 MyType *const *const U, MyType *const b, MyType *const x, const MyType omega,
                 const int n, const MyType *const exact)  // Метод релаксации
{
    MyType **B,
            **invB,
            **C,
            **C_L,
            **C_D,
            **C_U,
            *y,
            *x_k,
            *diff;

    B = new MyType *[n];
    invB = new MyType *[n];
    C = new MyType *[n];
    C_L = new MyType *[n];
    C_D = new MyType *[n];
    C_U = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        B[i] = new MyType[n];
        invB[i] = new MyType[n];
        C[i] = new MyType[n];
        C_L[i] = new MyType[n];
        C_D[i] = new MyType[n];
        C_U[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = D[i][j] + omega * L[i][j];
        }
    }

    int result = computeInverseMatrix(B, invB, n);
    if (result != 0) {
        cout << "Матрица B вырождена\n";
        return 1;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += invB[i][k] * A[k][j];
            }
            MyType e_ij = (i == j) ? 1 : 0;
            C[i][j] = e_ij - omega * sum;
        }
    }

    cout << "Матрица C: " << '\n';
    printMatrix(C, n);
    MyType cNorm = normMatrix(C, n);
    cout << "Норма матрицы С: " << cNorm << '\n';

    LDU(C, C_L, C_D, C_U, n);
    MyType clNorm = normMatrix(C_L, n);
    cout << "||C_L||_inf = " << clNorm << '\n';
    MyType cdNorm = normMatrix(C_D, n);
    cout << "||C_D||_inf = " << cdNorm << '\n';
    MyType cuNorm = normMatrix(C_U, n);
    cout << "||C_U||_inf = " << cuNorm << '\n';
    cout << "||C_L||_inf + ||C_U||_inf = " << clNorm + cuNorm << '\n';

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = (i == j) ? 0 : -A[i][j] / A[i][i];
        }
        y[i] = omega * b[i] / A[i][i];
        x_k[i] = y[i];
    }

    double iter_est = 0.0;
    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += C[i][j] * x_k[j];
            }
            x[i] = (1 - omega) * x_k[i] + omega * sum + y[i];
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += C[i][j] * x[j];
            }
            x[i] += omega * sum;
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
        if (iter == 1) {
            iter_est = log(EPS / normVector(diff, n) * (1 - cNorm)) / log(cNorm);
        }
    } while (normVector(diff, n) > (EPS * (1 - cNorm)) / cuNorm);

    cout << "Вектор y: ";
    printVector(y, n);
    cout << "Оценка числа итераций: > " << iter_est << '\n';
    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    printError(exact, x, n, normVector);
    printResidual(A, b, x, n, normVector);

    for (int i = 0; i < n; ++i) {
        delete[] B[i];
        delete[] invB[i];
        delete[] C[i];
        delete[] C_L[i];
        delete[] C_D[i];
        delete[] C_U[i];
    }
    delete[] B;
    delete[] invB;
    delete[] C;
    delete[] C_L;
    delete[] C_D;
    delete[] C_U;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int seidelThreeDiag(const MyType *const a, const MyType *const b, const MyType *const c,
                    const MyType *const d, MyType *const x, const MyType omega,
                    const int n)  // Метод релаксации для трехдиагональных матриц
{
    MyType *y,
            *x_k,
            *diff;

    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        y[i] = omega * d[i] / b[i];
        x_k[i] = y[i];
    }

    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType s = i < n ? c[i] / b[i] * x_k[i + 1] : 0.0;
            x[i] = (1 - omega) * x_k[i] - omega * s + y[i];
        }

        for (int i = 0; i < n; ++i) {
            MyType s = i > 0 ? a[i] / b[i] * x[i - 1] : 0.0;
            x[i] -= omega * s;
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
    } while (normVector(diff, n) > EPS);

    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);

    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}