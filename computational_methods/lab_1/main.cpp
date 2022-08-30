#include <cmath>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;

typedef double MyType;
const MyType eps0 = (typeid(MyType).name()[0] == 'd') ? 1e-12 : 1e-6;
const char *const settingsFileName = "settings.dat";

int readSettings(const char *fileName, size_t &n, string &matrixName, string &vectorName,
                 string &outputName);

void memoryAllocation(MyType **&A, MyType **&Q, MyType **&R, MyType **&T, MyType **&invA,
                      MyType *&b, MyType *&tb, MyType *&x, int n);

void memoryDeallocation(MyType **&A, MyType **&Q, MyType **&R, MyType **&T, MyType **&invA,
                        MyType *&b, MyType *&tb, MyType *&x, int n);

int readData(MyType *const *A, MyType *b, const string &matrixName, const string &vectorName);

void printSystem(const MyType *const *A, const MyType *b, int n);

void printSolution(const MyType *x, int n);

void printMatrix(const MyType *const *A, int n);

int gaussMethod(MyType **A, MyType *b, MyType *x, int n);

void backward(const MyType *const *A, const MyType *b, MyType *x, int n);

void printResidual(const MyType *const *A, const MyType *b, const MyType *x, int n,
                   MyType normVector(const MyType *const, const int));

int qrMethod(const MyType *const *A, const MyType *b, MyType *const *Q, MyType *const *R,
             MyType *const *T, MyType *tb, MyType *x, int n);

MyType computeNorm2(const MyType *x, int n);

MyType computeNorm1(const MyType *x, int n);

MyType computeNormInf(const MyType *x, int n);

MyType computeMatrixNorm1(const MyType *const *A, int n);

MyType computeMatrixNormInf(const MyType *const *A, int n);

void computeInverseMatrix(const MyType *const *R, const MyType *const *T, MyType *const *invA, int n);

MyType **multiplyMatrix(const MyType **A, const MyType **B, int n);

void conditionality(MyType *const *A, MyType *const *invA, MyType *const *R, MyType *const *T,
                    const MyType *b, MyType *x, int n,
                    MyType normVector(const MyType *const, const int),
                    MyType normMatrix(const MyType *const *const, const int),
                    MyType addNoise1(MyType), MyType addNoise2(MyType));

MyType addNoisePlus(MyType x);

MyType addNoiseMinus(MyType x);


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

    MyType **A,      // Матрица коэффициентов
            **Q,     // Ортоганальная матрица
            **R,     // Верхнетреугольная матрица
            **T,     // Матрица результирующего вращения
            **invA,  // Единичная матрица
            *b,      // Вектор правой части
            *tb,     // Вектор правой части, умноженный слева на матрицу вращения
            *x;      // Вектор неизвестных
    memoryAllocation(A, Q, R, T, invA, b, tb, x, n);

    result = readData(A, b, matrixFileName, vectorFileName);
    if (result != 0) {
        return result;
    }

    cout << "\nВариант 19, система " << matrixFileName[11] << ", тип "
         << ((typeid(MyType).name()[0] == 'd') ? "double" : "float") << "\n\n";
    printSystem((const MyType **) A, b, n);

    result = gaussMethod(A, b, x, n);
    if (result == 0) {
        printSolution(x, n);
        printResidual((const MyType **) A, b, x, n, computeNorm1);
        printResidual((const MyType **) A, b, x, n, computeNormInf);
    }

    result = readData(A, b, matrixFileName, vectorFileName);
    if (result != 0) {
        return result;
    }

    result = qrMethod((const MyType **) A, b, Q, R, T, tb, x, n);
    if (result == 0) {
        printSolution(x, n);
        printResidual((const MyType **) A, b, x, n, computeNorm1);
        printResidual((const MyType **) A, b, x, n, computeNormInf);

        cout << "\nМатрица Q:\n";
        printMatrix((const MyType **) Q, n);
        cout << "Матрица R:\n";
        printMatrix((const MyType **) R, n);

        cout << "\nВариант 19, система " << matrixFileName[11] << ", тип "
             << ((typeid(MyType).name()[0] == 'd') ? "double" : "float") << '\n';

        cout << "\nМатрица A^-1:\n";
        computeInverseMatrix((const MyType **) R, (const MyType **) T, invA, n);
        printMatrix((const MyType **) invA, n);
        cout << "Матрица A^-1 * A:\n";
        MyType **A1A = multiplyMatrix((const MyType **) invA, (const MyType **) A, n);
        printMatrix((const MyType **) A1A, n);
        for (int i = 0; i < n; ++i) {
            delete[] A1A[i];
        }
        delete[] A1A;

        cout << "\nОктаэдрическая норма:\n";
        conditionality(A, invA, R, T, b, x, n, computeNorm1, computeMatrixNorm1, addNoisePlus,
                       addNoiseMinus);
        cout << "\nКубическая норма:\n";
        conditionality(A, invA, R, T, b, x, n, computeNormInf, computeMatrixNormInf, addNoisePlus,
                       addNoiseMinus);
    }

    memoryDeallocation(A, Q, R, T, invA, b, tb, x, n);

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


void memoryAllocation(MyType **&A, MyType **&Q, MyType **&R, MyType **&T, MyType **&invA,
                      MyType *&b, MyType *&tb, MyType *&x, const int n)  // Выделение памяти
{
    A = new MyType *[n];
    Q = new MyType *[n];
    R = new MyType *[n];
    T = new MyType *[n];
    invA = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new MyType[n];
        Q[i] = new MyType[n];
        R[i] = new MyType[n];
        T[i] = new MyType[n];
        invA[i] = new MyType[n];
    }
    b = new MyType[n];
    tb = new MyType[n];
    x = new MyType[n];
}


void memoryDeallocation(MyType **&A, MyType **&Q, MyType **&R, MyType **&T, MyType **&invA,
                        MyType *&b, MyType *&tb, MyType *&x, const int n)  // Освобождение памяти
{
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] Q[i];
        delete[] R[i];
        delete[] T[i];
        delete[] invA[i];
    }
    delete[] A;
    delete[] Q;
    delete[] R;
    delete[] T;
    delete[] invA;
    delete[] b;
    delete[] tb;
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


void printSolution(const MyType *const x, const int n)  // Вывод вектора решения
{
    cout << "Вектор решения: { " << x[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << x[i];
    }
    cout << " }\n";
}


void printMatrix(const MyType *const *const A, const int n)  // Вывод матрицы
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << std::setw(15) << A[i][j];
        }
        cout << '\n';
    }
}


int gaussMethod(MyType **const A, MyType *const b, MyType *const x, const int n)  // Метода Гаусса
{
    cout << "\nМетод Гаусса:\n";

    MyType *temp;  // Вектор-буфер
    MyType t;      // Значение-буфер

    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n; ++k) {
        // Частичный поиск ведущего элемента по столбцу
        int iMax = k;
        for (int i = k; i < n; ++i) {
            if (fabs(A[i][k]) > fabs(A[iMax][k])) {
                iMax = k;
            }
        }

        // Если на диагонали нулевой элемент
        if (fabs(A[iMax][k]) < eps0) {
            cout << "Матрица вырождена. Система несовместна.\n"
                 << endl;
            return 1;
        }

        // Меняем строки
        if (iMax != k) {
            temp = A[k];
            A[k] = A[iMax];
            A[iMax] = temp;

            t = b[k];
            b[k] = b[iMax];
            b[iMax] = t;
        }

        // Прямой ход
        for (int i = k + 1; i < n; ++i) {
            t = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= A[k][j] * t;
            }
            b[i] -= b[k] * t;
        }
    }

    backward((const MyType **) A, b, x, n);

    return 0;
}


void backward(const MyType *const *const A, const MyType *const b, MyType *const x, const int n)  // Обратный ход
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
        if (fabs(x[i]) < eps0) {
            x[i] = 0.0;
        }
    }
}


void printResidual(const MyType *const *const A, const MyType *const b, const MyType *const x, const int n,
                   MyType normVector(const MyType *const, const int))  // Подсчет невязки
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


int qrMethod(const MyType *const *const A, const MyType *const b, MyType *const *const Q, MyType *const *const R,
             MyType *const *const T, MyType *const tb, MyType *const x, const int n)  // Метод QR-разложения
{
    cout << "\nМетод QR-разложения:\n";

    // Инициализируем матрицы
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            T[i][j] = (i == j) ? 1.0 : 0.0;
            R[i][j] = A[i][j];  // R = T*A
        }
    }

    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            if (R[i][k] == 0)
                continue;  // Вращение не требуется

            MyType denominator = sqrt(R[k][k] * R[k][k] + R[i][k] * R[i][k]);

            MyType c = R[k][k] / denominator;
            MyType s = R[i][k] / denominator;
            R[k][k] = denominator;  // Диагональный элемент будет равен знаменателю
            R[i][k] = 0.0;          // А этот должен обнулиться

            // Домножаем матрицу A слева на матрицу T_ij
            for (int j = k + 1; j < n; ++j) {
                MyType Rkj = R[k][j];
                MyType Rij = R[i][j];
                R[k][j] = c * Rkj + s * Rij;
                R[i][j] = -s * Rkj + c * Rij;
            }

            // Домножаем матрицу T слева на матрицу T_ij
            for (int j = 0; j < n; ++j) {
                MyType Tkj = T[k][j];
                MyType Tij = T[i][j];
                T[k][j] = c * Tkj + s * Tij;
                T[i][j] = -s * Tkj + c * Tij;
            }
        }
    }

    // Если на диагонали есть нулевой элемент
    for (int k = 0; k < n - 1; ++k) {
        if (R[k][k] < eps0) {
            cout << "Матрица вырождена. Система несовместна.\n"
                 << endl;
            return 1;
        }
    }

    // Домножаем вектор правой части b слева на матрицу T
    for (int i = 0; i < n; ++i) {
        MyType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += T[i][j] * b[j];
        }
        tb[i] = sum;
    }

    // Вычисляем матрицу Q транспонированием матрицы T
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q[i][j] = T[j][i];
        }
    }

    backward((const MyType **) R, tb, x, n);

    return 0;
}


MyType computeNorm2(const MyType *const x, const int n)  // Евклидова норма
{
    MyType norm = 0;
    for (int i = 0; i < n; ++i) {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
}


MyType computeNorm1(const MyType *const x, const int n)  // Октаэдрическая норма вектора
{
    MyType sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += fabs(x[i]);
    }
    return sum;
}


MyType computeNormInf(const MyType *const x, const int n)  // Кубическая норма вектора
{
    MyType max = 0.0;
    for (int i = 0; i < n; ++i) {
        if (fabs(x[i]) > max) {
            max = fabs(x[i]);
        }
    }
    return max;
}


MyType computeMatrixNorm1(const MyType *const *const A, const int n)  // Октаэдрическая норма матрицы
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


MyType computeMatrixNormInf(const MyType *const *const A, const int n)  // Кубическая норма матрицы
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


void computeInverseMatrix(const MyType *const *const R, const MyType *const *const T,
                          MyType *const *const invA, const int n)  // Вычисление обратой матрицы
{
    MyType **E = new MyType *[n];  // Единичная матрица
    for (int i = 0; i < n; ++i) {
        E[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            E[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    MyType *x_j = new MyType[n];   // j-й стобец искомой матрицы
    MyType *te_j = new MyType[n];  // j-й столбец матрицы E, умноженный слева на матрицу T

    // Решаем систему R * x_j = T * e_j
    for (int j = 0; j < n; ++j) {
        // Домножаем e_j слева на матрицу T
        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += T[i][k] * E[k][j];
            }
            te_j[i] = sum;
        }
        backward((const MyType **) R, te_j, x_j, n);

        for (int i = 0; i < n; ++i) {
            invA[i][j] = x_j[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        delete[] E[i];
    }
    delete[] E;
    delete[] x_j;
    delete[] te_j;
}


MyType **multiplyMatrix(const MyType **A, const MyType **B, int n)  // Умножение квадратных матриц
{
    MyType **M = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        M[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += A[i][k] * B[k][j];
            }
            M[i][j] = sum;
        }
    }
    return M;
}


void conditionality(MyType *const *const A, MyType *const *const invA, MyType *const *const R,
                    MyType *const *const T, const MyType *const b, MyType *const x, int n,
                    MyType normVector(const MyType *const, const int),
                    MyType normMatrix(const MyType *const *const, const int),
                    MyType addNoise1(MyType), MyType addNoise2(MyType))  // Число обусловленности и его оценка
{
    MyType *b1, *b2,         // Векторы правой части после возмущения
            *tb1, *tb2,      // Векторы правой части после возмущения, умноженные слева на матрицу T элементарных поворотов
            *x1, *x2,        // Векторы неизвестных после возмущения
            db, dx,          // Относительные погрешности b и x
            *Db, *Dx,        // Абсолютные погрешности b и x
            *cond1, *cond2,  // Числа обусловленности при раличных возмущениях
            condEstimate;    // Оценка обусловленности матрицы

    b1 = new MyType[n];
    b2 = new MyType[n];
    tb1 = new MyType[n];
    tb2 = new MyType[n];
    x1 = new MyType[n];
    x2 = new MyType[n];
    Db = new MyType[n];
    Dx = new MyType[n];
    cond1 = new MyType[n];
    cond2 = new MyType[n];

    // Для всех элементов вектора правой части:
    for (int k = 0; k < n; ++k) {
        // Добавляем возмущение
        b1[k] = addNoise1(b[k]);
        b2[k] = addNoise2(b[k]);

        // Домножаем слева на T
        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += T[i][j] * b1[j];
            }
            tb1[i] = sum;

            sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += T[i][j] * b2[j];
            }
            tb2[i] = sum;
        }

        // Получем решение соответствующей системы
        backward((const MyType **) R, tb1, x1, n);
        backward((const MyType **) R, tb2, x2, n);

        // Считаем оценку числа обусловленности
        for (int i = 0; i < n; ++i) {
            Dx[i] = x1[i] - x[i];
            Db[i] = b1[i] - b[i];
        }
        dx = normVector(Dx, n) / normVector(x, n);
        db = normVector(Db, n) / normVector(b, n);
        cond1[k] = dx / db;

        for (int i = 0; i < n; ++i) {
            Dx[i] = x2[i] - x[i];
            Db[i] = b2[i] - b[i];
        }
        dx = normVector(Dx, n) / normVector(x, n);
        db = normVector(Db, n) / normVector(b, n);
        cond2[k] = dx / db;

        // Возвращаем значение до возмущения
        b1[k] = b[k];
        b2[k] = b[k];
    }

    MyType max1 = 0.0, max2 = 0.0;
    int iMax1, iMax2, iMax;  // Номера максимумов

    // Выбираем максимальную оценку
    for (int i = 0; i < n; ++i) {
        if (cond1[i] > max1) {
            max1 = cond1[i];
            iMax1 = i;
        }
        if (cond2[i] > max2) {
            max2 = cond2[i];
            iMax2 = i;
        }
    }

    if (max1 >= max2) {
        condEstimate = max1;
        iMax = iMax1;
        b1[iMax] = addNoise1(b[iMax]);

        for (int k = 0; k < n; ++k) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += T[k][j] * b1[j];
            }
            tb1[k] = t;
        }

        backward((const MyType **) R, tb1, x1, n);
        printResidual((const MyType **) R, tb1, x1, n, normVector);
    } else {
        condEstimate = max2;
        iMax = iMax2;
        b2[iMax] = addNoise2(b[iMax]);

        for (int k = 0; k < n; ++k) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += T[k][j] * b2[j];
            }
            tb2[k] = t;
        }

        backward((const MyType **) R, tb2, x2, n);
        printResidual((const MyType **) R, tb2, x2, n, normVector);
    }

    cout << "Номер компоненты вектора правой части, наиболее сильно влияющей на решение: "
         << iMax << '\n';
    cout << "Оценка числа обусловленности: >= " << condEstimate << '\n';

    MyType normA = normMatrix((const MyType **) A, n);
    cout << "Норма матрицы: " << normA << '\n';
    MyType normInvA = normMatrix((const MyType **) invA, n);
    cout << "Норма обратной матрицы: " << normInvA << '\n';
    cout << "Обусловленность матрицы: " << normA * normInvA << '\n';

    delete[] tb1;
    delete[] tb2;
    delete[] b1;
    delete[] b2;
    delete[] x1;
    delete[] x2;
    delete[] Dx;
    delete[] Db;
    delete[] cond1;
    delete[] cond2;
}


MyType addNoisePlus(const MyType x)  // Добавляем 0.01
{
    return x + 0.01;
}


MyType addNoiseMinus(const MyType x)  // Вычитаем 0.01
{
    return x - 0.01;
}
