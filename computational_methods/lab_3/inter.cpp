#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>

int NumIntNodes;  // Количество узлов интерполяции
const double PI = M_PI;
double eps = 0.01;  // 0.0001


// Интерполируемая функция
double Func(double x) {
    // return exp(x);
    // return x * x;
    // return 1 / (1 + x * x);
    // return 1 / atan(1 + 10 * x * x);// [-3,3]
    // return pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
    // return 1;
    // return x;
    // return 1 / (1 + 25 * x * x);
    // return sin(PI * x);// [-1,1], [-1.25,1.25]
    if (x < 0)  // [-1,1], [-1.25,1.25]
        return PI * x;
    else
        return sin(PI * x);
}

// Печать вектора
void PrintVector(double* vec) {
    if (vec != nullptr) {
        for (auto i = 0; i < NumIntNodes; ++i) {
            std::cout << std::setw(7) << vec[i] << "\n";
        }
        std::cout << "\n";
    } else {
        std::cout << "Print error: nullptr!\n";
    }
}

// Запись значений интерполянта
void PrintToFile(const std::string& filepath, int size, double* xValues, double* fValues) {
    std::ofstream fout;
    fout.open(filepath);
    if (!fout.is_open()) {
        std::cout << "Error file!\n";
    } else {
        for (int i = 0; i <= size; i++) {
            fout << xValues[i] << " " << fValues[i] << "\n";
        }
    }
    fout.close();
    fout.clear();
}

// Построение равномерной сетки
void UniformGrid(int NumNodes, double* xValues, double a, double b) {
    double h = (b - a) / NumNodes;
    for (auto i = 0; i <= NumNodes; ++i)
        xValues[i] = a + h * i;
}

// Построение чебышевской сетки
void ChebyshevGrid(double* xValues, double a, double b) {
    for (auto i = 0; i <= NumIntNodes; ++i)
        xValues[NumIntNodes - i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * PI / (2 * (NumIntNodes + 1)));
}

// Вычисление значений функции в узлах сетки
void FuncValues(int NumNodes, double* xValues, double* fValues) {
    for (auto i = 0; i <= NumNodes; ++i)
        fValues[i] = Func(xValues[i]);
}

// Норма ошибки интерполяции
double ErrorNorm(int InterpSize, double* xinterpValues, double* finterpValues) {
    double* fValues_err = new double[InterpSize + 1];
    FuncValues(InterpSize, xinterpValues, fValues_err);
    double errNorm = 0;
    for (auto i = 0; i <= InterpSize; ++i) {
        if (fabs(fValues_err[i] - finterpValues[i]) > errNorm) {
            errNorm = fabs(fValues_err[i] - finterpValues[i]);
        }
    }
    return errNorm;
}

// Правая 3-х диагональная прогонка
void RightTridiagRun(double* x, const double* a, const double* b, const double* c, const double* d) {
    double* alpha = new double[NumIntNodes];
    double* beta = new double[NumIntNodes];
    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (auto i = 1; i < NumIntNodes; ++i) {
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        beta[i] = (d[i] + a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }
    x[NumIntNodes - 1] = (d[NumIntNodes - 1] + a[NumIntNodes - 1] * beta[NumIntNodes - 2]) / (b[NumIntNodes - 1] - a[NumIntNodes - 1] * alpha[NumIntNodes - 2]);
    for (auto i = NumIntNodes - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    delete[] alpha;
    delete[] beta;
}

// Левая 3-х диагональная прогонка
void LeftTridiagRun(double* x, const double* a, const double* b, const double* c, const double* d) {
    double* xi = new double[NumIntNodes];
    double* eta = new double[NumIntNodes];
    xi[NumIntNodes - 1] = a[NumIntNodes - 1] / b[NumIntNodes - 1];
    eta[NumIntNodes - 1] = d[NumIntNodes - 1] / b[NumIntNodes - 1];
    for (auto i = NumIntNodes - 2; i > 0; --i) {
        xi[i] = a[i] / (b[i] - c[i] * xi[i + 1]);
        eta[i] = (d[i] + c[i] * eta[i + 1]) / (b[i] - c[i] * xi[i + 1]);
    }
    x[0] = 0;  // (d[0] + c[0] * eta[1]) / (b[0] - c[0] * xi[1])
    for (auto i = 1; i < NumIntNodes; ++i) {
        x[i] = xi[i] * x[i - 1] + eta[i];
    }
    x[NumIntNodes] = 0;
    delete[] xi;
    delete[] eta;
}

// Вычисление коэффициентов a,b,d
void CubicPolyCoeff(const double* xValues, const double* fValues, double* a, double* b, double* c, double* d) {
    double* h = new double[NumIntNodes];
    double* g = new double[NumIntNodes];

    for (auto i = 0; i < NumIntNodes; ++i) {
        a[i] = fValues[i];
        h[i] = xValues[i + 1] - xValues[i];
        g[i] = (fValues[i + 1] - fValues[i]) / h[i];
    }
    // Коэффициенты СЛАУ в прогонке
    double* A = new double[NumIntNodes];
    double* B = new double[NumIntNodes];
    double* C = new double[NumIntNodes];
    double* D = new double[NumIntNodes];
    for (auto i = 1; i < NumIntNodes; ++i) {
        A[i] = h[i - 1];
        B[i] = -2 * (h[i - 1] + h[i]);
        C[i] = h[i];
        D[i] = -3 * (g[i] - g[i - 1]);
    }
    LeftTridiagRun(c, A, B, C, D);
    for (auto i = 0; i < NumIntNodes; ++i) {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] h;
    delete[] g;
}

// Сплайн-интерполяция
void SplineInterpolation(double* xValues, double* fValues, double leftb, double rightb) {
    double* a = new double[NumIntNodes];
    double* b = new double[NumIntNodes];
    double* c = new double[NumIntNodes + 1];
    double* d = new double[NumIntNodes];
    CubicPolyCoeff(xValues, fValues, a, b, c, d);

    // Значения интерполянта на мелкой сетке
    int InterpSize = NumIntNodes * 10;
    double* xinterpValues = new double[InterpSize + 1];
    double* finterpValues = new double[InterpSize + 1];
    UniformGrid(InterpSize, xinterpValues, leftb, rightb);
    int j = 0;
    double tempDifference;
    for (auto i = 0; i < NumIntNodes; ++i) {
        while ((xinterpValues[j] <= xValues[i + 1]) && (j <= InterpSize)) {
            tempDifference = xinterpValues[j] - xValues[i];
            finterpValues[j] = a[i] + b[i] * tempDifference + c[i] * tempDifference * tempDifference + d[i] * tempDifference * tempDifference * tempDifference;
            ++j;
        }
    }
    while ((xinterpValues[j] > xValues[NumIntNodes]) && (j <= InterpSize)) {
        tempDifference = xinterpValues[j] - xValues[NumIntNodes - 1];
        finterpValues[j] = a[NumIntNodes - 1] + b[NumIntNodes - 1] * tempDifference + c[NumIntNodes - 1] * tempDifference * tempDifference + d[NumIntNodes - 1] * tempDifference * tempDifference * tempDifference;
        ++j;
    }

    PrintToFile("SplineInterpolation.dat", InterpSize, xinterpValues, finterpValues);
    double errNorm = ErrorNorm(InterpSize, xinterpValues, finterpValues);
    std::cout << "Error norm (spline interpolation) = " << errNorm << "\n";

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] xinterpValues;
    delete[] finterpValues;
}

// Интерполяционный полином Лагранжа в точке x
double LagrangePolynomial(double x, const double* xValues, const double* fValues) {
    double LagrPol = 0;
    double basePol;

    for (int k = 0; k <= NumIntNodes; k++) {
        basePol = 1;
        for (int i = 0; i <= NumIntNodes; i++) {
            if (i != k) {
                basePol *= (x - xValues[i]) / (xValues[k] - xValues[i]);
            }
        }
        LagrPol += basePol * fValues[k];
    }
    return LagrPol;
}

// Интерполяция полиномом Лагранджа
void LagrangianInterpolation(double* xValues, double* fValues, double leftb, double rightb) {
    // Значения интерполянта на мелкой сетке
    int InterpSize = NumIntNodes * 10;
    double* xinterpValues = new double[InterpSize + 1];
    double* finterpValues = new double[InterpSize + 1];
    UniformGrid(InterpSize, xinterpValues, leftb, rightb);
    for (auto i = 0; i <= InterpSize; ++i) {
        finterpValues[i] = LagrangePolynomial(xinterpValues[i], xValues, fValues);
    }

    PrintToFile("LagrangianInterpolation.dat", InterpSize, xinterpValues, finterpValues);
    double errNorm = ErrorNorm(InterpSize, xinterpValues, finterpValues);
    std::cout << "Error norm (Lagrangian interpolation) = " << errNorm << "\n";

    delete[] xinterpValues;
    delete[] finterpValues;
}

int main() {
    char stopCh;
    char Grid;             // Выбор сетки
    double leftb, rightb;  // Интервал интерполирования

    do {
        std::cout << "Number of interpolation nodes = ";
        std::cin >> NumIntNodes;
        std::cout << "Interpolation interval: ";
        std::cin >> leftb >> rightb;
        std::cout << "Grid: a) Uniform grid, b) Chebyshev grid: ";
        std::cin >> Grid;
        double* xValues = new double[NumIntNodes + 1];
        double* fValues = new double[NumIntNodes + 1];
        switch (Grid) {
            case 'a':
                UniformGrid(NumIntNodes, xValues, leftb, rightb);
                break;
            case 'b':
                ChebyshevGrid(xValues, leftb, rightb);
                break;
            default:
                std::cout << "Incorrect grid type!";
                break;
        }
        FuncValues(NumIntNodes, xValues, fValues);
        SplineInterpolation(xValues, fValues, leftb, rightb);
        LagrangianInterpolation(xValues, fValues, leftb, rightb);
        // std::cout<<LagrangePolynomial(2.2, xValues, fValues);
        std::cout << "Continue? ";
        std::cin >> stopCh;
    } while (stopCh != '0');

    system("PAUSE");
    return 0;
}
