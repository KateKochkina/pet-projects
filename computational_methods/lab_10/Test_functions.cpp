#include "Test_functions.h"

//test 1
double K1(double x, double s) {
	return (1.0 - x * cos(x * s)) / 2.0;
}

double f1(double x) {
	return (1.0 + sin(x)) / 2.0;
}

double g1(double x) {
	return x * x + sqrt(x);
}

//Taylor decomposition
double phi0(double x) {
	return 1.0 / 2;
}
double psi0(double s) {
	return 1.0;
}

double phi1(double x) {
	return -x / 2;
}
double psi1(double s) {
	return 1.0;
}

double phi2(double x) {
	return pow(x, 3) / 4;
}
double psi2(double s) {
	return pow(s, 2);
}

double phi3(double x) {
	return -pow(x, 5) / 24;
}
double psi3(double s) {
	return pow(s, 4);
}

double phi4(double x) {
	return pow(x, 7) / 1440;
}
double psi4(double s) {
	return pow(s, 6);
}

double phi5(double x) {
	return -pow(x, 9) / 80640;
}
double psi5(double s) {
	return pow(s, 8);
}

double phi6(double x) {
	return pow(x, 11) / 7257600;
}
double psi6(double s) {
	return pow(s, 10);
}

double phi7(double x) {
	return -pow(x,13) / 958003200;
}
double psi7(double s) {
	return pow(s,12);
}

double phi8(double x) {
	return pow(x, 15) / 174356582400;
}
double psi8(double s) {
	return pow(s, 14);
}

//17 var
double K17(double x, double s) {
	return x * (exp(-x * s) - 1.0);
}

double f17(double x) {
	return x + exp(-x);
}

double g17(double x) {
	return log2(x + 1.0) + x / 2.0;
}

//Taylor decomposition
double phi0_var17(double x) {
    return -pow(x, 2);
}
double psi0_var17(double s) {
    return s;
}

double phi1_var17(double x) {
    return pow(x, 3) / 2;
}
double psi1_var17(double s) {
    return pow(s, 2);
}

double phi2_var17(double x) {
    return -pow(x, 4) / 6;
}
double psi2_var17(double s) {
    return pow(s, 3);
}

double phi3_var17(double x) {
    return pow(x, 5) / 24;
}
double psi3_var17(double s) {
    return pow(s, 4);
}

double phi4_var17(double x) {
    return -pow(x, 6) / 120;
}
double psi4_var17(double s) {
    return pow(s, 5);
}

double phi5_var17(double x) {
    return pow(x, 7) / 720;
}
double psi5_var17(double s) {
    return pow(s, 6);
}

double phi6_var17(double x) {
    return -pow(x, 8) / 5040;
}
double psi6_var17(double s) {
    return pow(s, 7);
}

double phi7_var17(double x) {
    return pow(x, 9) / 40320;
}
double psi7_var17(double s) {
    return pow(s, 8);
}

double phi8_var17(double x) {
    return -pow(x, 10) / 362880;
}
double psi8_var17(double s) {
    return pow(s, 9);
}

//singular equation
TVector<double> Q_sing(TVector<double> r, TVector<double> rho) {
	double znam = 2 * PI * ((r[0] - rho[0]) * (r[0] - rho[0]) + (r[1] - rho[1]) * (r[1] - rho[1]));
	return { -(r[1] - rho[1]) / znam, (r[0] - rho[0]) / znam };
}

//s = phi * r => phi = s /r
double f_odd(int N, double phi) {
	return sin((N + 1) * phi / 2);
}

double f_even(int N, double phi) {
	return cos(N * phi / 2);
}
