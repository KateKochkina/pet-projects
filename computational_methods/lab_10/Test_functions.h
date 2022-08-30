#pragma once
#include<iostream>
#include"TVector.h"

const double PI = 3.14159265358979323846;


//test 1
double K1(double x, double s);

double f1(double x);

double g1(double x);

//Taylor decomposition
double phi0(double x);
double psi0(double s);
double phi1(double x);
double psi1(double s);
double phi2(double x);
double psi2(double s);
double phi3(double x);
double psi3(double s);
double phi4(double x);
double psi4(double s);
double phi5(double x);
double psi5(double s);
double phi6(double x);
double psi6(double s);
double phi7(double x);
double psi7(double s);
double phi8(double x);
double psi8(double s);

//17 var
double K17(double x, double s);

double f17(double x);

double g17(double x);

double phi0_var17(double x);
double psi0_var17(double s);
double phi1_var17(double x);
double psi1_var17(double s);
double phi2_var17(double x);
double psi2_var17(double s);
double phi3_var17(double x);
double psi3_var17(double s);
double phi4_var17(double x);
double psi4_var17(double s);
double phi5_var17(double x);
double psi5_var17(double s);
double phi6_var17(double x);
double psi6_var17(double s);
double phi7_var17(double x);
double psi7_var17(double s);
double phi8_var17(double x);
double psi8_var17(double s);

//singular equation
TVector<double> Q_sing(TVector<double> r, TVector<double> rho);

double f_odd(int N, double phi);

double f_even(int N, double phi);
