#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TVector.h"
#include "TMatrix.h"

const double eps = 1e-16; //точность для метода Гаусса

template <typename MYTYPE> void filing_gauss(std::string filename, TVector<MYTYPE> vec) {
	std::ofstream fout;
	fout.open(filename);
	fout << vec.get_dim() << '\n';
	for (auto j = 0; j < vec.get_dim(); ++j) {
		fout << std::setprecision(20) << std::fixed << vec[j] << '\n';
	}
	fout << std::endl;
	fout.close();
	fout.clear();
}

void clear_file_gauss(std::string filename);

//Печать для метода Гаусса
template <typename MYTYPE> void PrintGauss(TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<int>& NRows, TVector<int>& NCols, int dim)
{
	std::cout << "A: \n";
	for (auto i = 0; i < dim; i++) {
		for (auto j = 0; j < dim; j++) {
			std::cout << std::setw(12) << A[NRows[i] * dim + NCols[j]] << " ";
		}
		std::cout << "\n";
	}

	std::cout << "\nb:\n";
	for (auto i = 0; i < dim; ++i) {
		std::cout << std::setw(7) << b[NRows[i]] << "\n";
	}
	std::cout << "\n";
}
//Проверка вырожденности
template <typename MYTYPE> bool DegeneracyMatr(MYTYPE MaxElem)
{
	if (std::fabs(MaxElem) < eps) return true;
	else return false;
}
//Выбор главного элемента
template <typename MYTYPE> bool MainElement(int iter, TMatrix<MYTYPE>& A, TVector<int>& NRows, TVector<int>& NCols, int dim)
{
	double MaxElem = std::fabs(A(NRows[iter], NCols[iter]));
	int MaxElemPos[2] = { iter,iter };
	bool transposition = false;

	//Поиск наибольшего по модулю элемента
	for (auto j = iter; j < dim; ++j) {
		for (auto k = iter; k < dim; ++k) {
			if (std::fabs(A(NRows[j], NCols[k])) > std::fabs(MaxElem)) {
				MaxElem = A(NRows[j], NCols[k]);
				MaxElemPos[0] = j;
				MaxElemPos[1] = k;
				transposition = true;
			}
		}
	}

	if (!DegeneracyMatr(MaxElem)) {
		//Перестановка строк/столбцов
		if (transposition == true) {
			int iRow = NRows[iter];
			int iCol = NCols[iter];
			NRows[iter] = NRows[MaxElemPos[0]];
			NCols[iter] = NCols[MaxElemPos[1]];
			NRows[MaxElemPos[0]] = iRow;
			NCols[MaxElemPos[1]] = iCol;
		}
		return true;
	}
	else return false;
}
//Прямой ход
template <typename MYTYPE> bool DirectGauss(TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<int>& NRows, TVector<int>& NCols)
{
	int dim = b.get_dim();
	for (auto iter = 0; iter < dim; ++iter) {
		if (MainElement(iter, A, NRows, NCols, dim)) {
			for (auto j = iter + 1; j < dim; ++j) {
				auto cft = A(NRows[iter], NCols[j]) / A(NRows[iter], NCols[iter]);
				//std::cout << "1\n";
				for (auto k = iter + 1; k < dim; ++k) {
					//std::cout << "1 + 2\n";
					A(NRows[k], NCols[j]) -= A(NRows[k], NCols[iter]) * cft;
				}
			}
			auto cft2 = b[NRows[iter]] / A(NRows[iter], NCols[iter]);
			//std::cout << "1\n";
			for (auto k = iter + 1; k < dim; ++k) {
				//std::cout << "1 + 2\n";
				b[NRows[k]] -= A(NRows[k], NCols[iter]) * cft2;
				A(NRows[k], NCols[iter]) = MYTYPE(0);
			}
		}
		else return false;
	}
	return true;
}
//Обратный ход (для вызова из GaussMeth)
template <typename MYTYPE> void ReverseGauss(TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x, TVector<int>& NRows, TVector<int>& NCols)
{
	int dim = b.get_dim();
	double sum;
	for (auto i = dim - 1; i >= 0; --i) {
		sum = 0;
		for (auto j = dim - 1; j >= i + 1; --j) {
			//std::cout << "1 + 2\n";
			sum += A(NRows[i], NCols[j]) * x[NCols[j]];
		}
		x[NCols[i]] = (b[NRows[i]] - sum) / A(NRows[i], NCols[i]);
		//std::cout << "1 + 2\n";
	}
}
//Метод Гаусса
template <typename MYTYPE> void GaussMeth(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, TVector<MYTYPE>& x, char print = 'n')
{
	bool degenerate = false;
	int dim; //размерность системы
	dim = b.get_dim();
	TVector<int> NRows(dim);
	TVector<int> NCols(dim);
	for (auto i = 0; i < dim; i++) {
		NRows[i] = i;
		NCols[i] = i;
	}
	//TMatrix<MYTYPE> A(matr);
	//TVector<MYTYPE> b(vec);
	//TVector<MYTYPE> x(dim);

	if (DirectGauss(A, b, NRows, NCols)) {
		ReverseGauss(A, b, x, NRows, NCols);
		if (print != 'n') {
			//std::cout << "Result of the Gauss method:\n";
			//PrintGauss(A, b, NRows, NCols);
			//std::cout << "x:\n" << x;
			//clear_file_gauss("realSOL_" + filename + ".txt");
			//filing_gauss("realSOL_" + filename + ".txt", x);
			clear_file_gauss(filename);
			filing_gauss(filename, x);
		}
	}
	else {
		degenerate = true;
		x.zero_vector();
		std::cout << "Degenerate matrix!\n";
	}
}

//Метод Гаусса (с возвратом)
template <typename MYTYPE> TVector<MYTYPE> GaussMeth(std::string filename, TMatrix<MYTYPE>& A, TVector<MYTYPE>& b, char print = 'n')
{
	bool degenerate = false;
	int dim; //размерность системы
	dim = b.get_dim();
	TVector<int> NRows(dim);
	TVector<int> NCols(dim);
	for (auto i = 0; i < dim; i++) {
		NRows[i] = i;
		NCols[i] = i;
	}
	TVector<MYTYPE> x(dim);

	if (DirectGauss(A, b, NRows, NCols)) {
		ReverseGauss(A, b, x, NRows, NCols);
		if (print != 'n') {
			clear_file_gauss(filename);
			filing_gauss(filename, x);
		}
	}
	else {
		degenerate = true;
		x.zero_vector();
		std::cout << "Degenerate matrix!\n";
	}
	return x;
}
