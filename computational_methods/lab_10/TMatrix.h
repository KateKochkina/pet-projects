#pragma once
#include "TVector.h"
template <typename MYTYPE> class TMatrix
{
protected:
	int m_nrows;
	int m_ncols;
	int m_elem_num;
	MYTYPE* m_pData;

	void init();
	void init(std::string filename);
	void close();
public:
	TMatrix() : m_nrows{ 0 }, m_ncols{ 0 }, m_elem_num{ 0 }, m_pData{ nullptr }{

	}
	TMatrix(const TMatrix<MYTYPE>& matr);
	TMatrix(int nrows, int ncols, MYTYPE* pData);
	TMatrix(int nrows, int ncols);
	TMatrix(std::string filename);
	TMatrix(TMatrix<MYTYPE>&& matr);
	~TMatrix() {
		close();
	}

	void copy(int nrows, int ncols, MYTYPE* pData = nullptr);
	void copy(const TMatrix<MYTYPE>& matr);
	void output();
	double norm(char norm_type);
	void add(const TMatrix<MYTYPE>& matr);
	void subtract(const TMatrix<MYTYPE>& matr);

	void zero_matrix();
	void identity_matrix();
	TMatrix<MYTYPE> transpose();
	//TMatrix<MYTYPE> precond_matr1();
	TMatrix<MYTYPE> precond_matr2();

	int get_nrows() const;
	int get_ncols() const;

	template <typename MYTYPE2> friend TMatrix<MYTYPE> matrix_product(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);
	template <typename MYTYPE2> friend TMatrix<MYTYPE> sum(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);
	template <typename MYTYPE2> friend TMatrix<MYTYPE> diff(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);

	TMatrix<MYTYPE>& operator += (const TMatrix<MYTYPE>& matr);
	TMatrix<MYTYPE>& operator -= (const TMatrix<MYTYPE>& matr);

	template <typename MYTYPE2> friend TMatrix<MYTYPE> operator + (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);
	template <typename MYTYPE2> friend TMatrix<MYTYPE> operator - (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);
	template <typename MYTYPE2> friend TMatrix<MYTYPE> operator * (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2);
	template <typename MYTYPE2> friend std::ostream& operator << (std::ostream& out, const TMatrix<MYTYPE>& matr);

	MYTYPE& operator () (const int row_ind, const int col_ind);
	const MYTYPE& operator () (const int row_ind, const int col_ind) const;

	MYTYPE& operator [] (const int index);
	const MYTYPE& operator [] (const int index) const;

	TMatrix<MYTYPE> operator - ();
	TMatrix<MYTYPE>& operator = (const TMatrix<MYTYPE>& matr);
	TMatrix<MYTYPE>& operator = (TMatrix<MYTYPE>&& matr);

	template <typename MYTYPE2> friend TVector<MYTYPE> operator * (const TMatrix<MYTYPE>& matr, const TVector<MYTYPE>& vec);
	template <typename MYTYPE2> friend TMatrix<MYTYPE> operator * (const MYTYPE numb, const TMatrix<MYTYPE>& matr);

	void set_as_column(const TVector<MYTYPE>& vec, int num_col);
	TVector<MYTYPE> get_column(int num_col);
	void remove_last_cols(int num_of_cols);
	void remove_last_rows(int num_of_rows);

	TMatrix(std::initializer_list<std::initializer_list<MYTYPE>> lst) : m_nrows(lst.size()), m_ncols(lst.begin()->size()), 
										m_elem_num(lst.size()* lst.begin()->size()), m_pData(new MYTYPE[lst.size() * lst.begin()->size()])
	{
		for (auto i = 0; i < m_nrows; ++i) {
			std::copy((lst.begin() + i)->begin(), (lst.begin() + i)->end(), m_pData + i * m_ncols);
		}
	}

};

template <typename MYTYPE> void TMatrix<MYTYPE>::init() {
	m_nrows = 0;
	m_ncols = 0;
	m_elem_num = 0;
	m_pData = nullptr;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::init(std::string filename) {
	std::ifstream fin;
	fin.open(filename, std::ios::in);
	if (!fin.is_open())
	{
		std::cout << "Error file!\n";
	}
	else
	{
		fin >> m_nrows >> m_ncols;
		m_elem_num = m_nrows * m_ncols;
		m_pData = new MYTYPE[m_elem_num];
		for (auto i = 0; i < m_elem_num; ++i) {
			fin >> m_pData[i];
		}
	}
	fin.close();
}

template <typename MYTYPE> void TMatrix<MYTYPE>::close() {
	if (m_nrows) delete[] m_pData;
	init();
}

template <typename MYTYPE> void TMatrix<MYTYPE>::copy(int nrows, int ncols, MYTYPE* pData) {
	close();
	if (nrows > 0) {
		m_nrows = nrows;
		m_ncols = ncols;
		m_elem_num = m_nrows * m_ncols;
		m_pData = new MYTYPE[m_elem_num];
		if (pData) {
			for (auto i = 0; i < m_elem_num; ++i) {
				m_pData[i] = pData[i];
			}
		}
	}
}

template <typename MYTYPE> void TMatrix<MYTYPE>::copy(const TMatrix<MYTYPE>& matr) {
	close();
	if (matr.m_nrows > 0) {
		m_nrows = matr.m_nrows;
		m_ncols = matr.m_ncols;
		m_elem_num = m_nrows * m_ncols;
		m_pData = new MYTYPE[m_elem_num];
		if (m_pData) {
			for (auto i = 0; i < m_elem_num; ++i) {
				m_pData[i] = matr.m_pData[i];
			}
		}
	}
}

template <typename MYTYPE> TMatrix<MYTYPE>::TMatrix(int nrows, int ncols, MYTYPE* pData) {
	copy(nrows, ncols, pData);
}

template <typename MYTYPE> TMatrix<MYTYPE>::TMatrix(int nrows, int ncols) {
	m_nrows = nrows;
	m_ncols = ncols;
	m_elem_num = m_nrows * m_ncols;
	m_pData = new MYTYPE[m_elem_num];
	for (auto i = 0; i < m_elem_num; ++i) {
		m_pData[i] = MYTYPE(0);
	}
}

template <typename MYTYPE> TMatrix<MYTYPE>::TMatrix(std::string filename) {
	init(filename);
}

template <typename MYTYPE> TMatrix<MYTYPE>::TMatrix(const TMatrix<MYTYPE>& matr) {
	copy(matr);
}

template <typename MYTYPE> TMatrix<MYTYPE>::TMatrix(TMatrix<MYTYPE>&& matr) {
	m_nrows = matr.m_nrows;
	m_ncols = matr.m_ncols;
	m_elem_num = matr.m_elem_num;
	m_pData = matr.m_pData;
	matr.m_nrows = 0;
	matr.m_ncols = 0;
	matr.m_elem_num = 0;
	matr.m_pData = nullptr;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::output() {
	for (auto i = 0; i < m_nrows; ++i) {
		for (auto j = 0; j < m_ncols; ++j) {
			std::cout << std::setw(12) << m_pData[i * m_ncols + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

template <typename MYTYPE> TMatrix<MYTYPE> matrix_product(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	TMatrix<MYTYPE> matr_res(matr1.m_nrows, matr2.m_ncols);
	if (matr1.m_ncols == matr2.m_nrows) {
		for (auto i = 0; i < matr_res.m_nrows; ++i) {
			for (auto j = 0; j < matr_res.m_ncols; ++j) {
				for (auto k = 0; k < matr1.m_ncols; ++k) {
					matr_res.m_pData[i * matr_res.m_ncols + j] += matr1.m_pData[i * matr1.m_ncols + k] * matr2.m_pData[k * matr2.m_ncols + j];
				}
			}
		}
	}
	else {
		std::cout << "Impossible to multiply matrices! The matrix product := zero matrix!\n";
	}
	return matr_res;
}

template <typename MYTYPE> double TMatrix<MYTYPE>::norm(char norm_type) {
	MYTYPE norm = MYTYPE(0);
	MYTYPE sum; //����� ��������� � ������� (��� '1') ��� � ������ (��� '3')
	switch (norm_type) 
	{
	case '1': //�������������� �����
		for (auto j = 0; j < m_ncols; ++j) {
			sum = 0;
			for (auto i = 0; i < m_nrows; ++i) {
				sum += std::abs(m_pData[i * m_ncols + j]);
			}
			if (sum > norm) {
				norm = sum;
			}
		}
		break;
	case '3': //���������� �����
		for (auto i = 0; i < m_nrows; ++i) {
			sum = 0;
			for (auto j = 0; j < m_ncols; ++j) {
				sum += std::abs(m_pData[i * m_ncols + j]);
			}
			if (sum > norm) {
				norm = sum;
			}
		}
		break;
	default: //�������� ��� �����
		std::cout << "Wrong norm type! Norm := 0!";
		break;
	}
	return(norm);
}

template <typename MYTYPE> TMatrix<MYTYPE> sum(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	TMatrix<MYTYPE> matr_res(matr1);
	if (matr1.m_nrows == matr2.m_nrows && matr1.m_ncols == matr2.m_ncols) {
		for (auto i = 0; i < matr_res.m_elem_num; ++i) {
			matr_res.m_pData[i] += matr2.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of matrix! The result is the first matrix!\n";
	}
	return matr_res;
}

template <typename MYTYPE> TMatrix<MYTYPE> diff(const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	TMatrix<MYTYPE> matr_res(matr1);
	if (matr1.m_nrows == matr2.m_nrows && matr1.m_ncols == matr2.m_ncols) {
		for (auto i = 0; i < matr_res.m_elem_num; ++i) {
			matr_res.m_pData[i] -= matr2.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of matrix! The result is the first matrix!\n";
	}
	return matr_res;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::add(const TMatrix<MYTYPE>& matr) {
	if (m_nrows == matr.m_nrows && m_ncols == matr.m_ncols) {
		for (auto i = 0; i < m_elem_num; ++i) {
			m_pData[i] += matr.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of matrtix! The original matrix is not changed!";
	}
}

template <typename MYTYPE> void TMatrix<MYTYPE>::subtract(const TMatrix<MYTYPE>& matr) {
	if (m_nrows == matr.m_nrows && m_ncols == matr.m_ncols) {
		for (auto i = 0; i < m_elem_num; ++i) {
			m_pData[i] -= matr.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of matrix! The original matrix is not changed!";
	}
}

template <typename MYTYPE> void TMatrix<MYTYPE>::zero_matrix() {
	for (auto i = 0; i < m_nrows * m_ncols; ++i) {
		m_pData[i] = MYTYPE(0);
	}
}

template <typename MYTYPE> void TMatrix<MYTYPE>::identity_matrix() {
	if (m_nrows == m_ncols) {
		for (auto i = 0; i < m_nrows; ++i) {
			for (auto j = 0; j < m_ncols; ++j) {
				if (i == j) {
					m_pData[i * m_ncols + j] = MYTYPE(1);
				}
				else {
					m_pData[i * m_ncols + j] = MYTYPE(0);
				}
			}
		}
	}
	else {
		std::cout << "It's not a square matrix! Impossible to make a single matrix!\n";
	}
}

template <typename MYTYPE> TMatrix<MYTYPE> TMatrix<MYTYPE>::transpose() {
	TMatrix<MYTYPE> matr_res(m_ncols, m_nrows);
	for (auto i = 0; i < matr_res.m_nrows; ++i) {
		for (auto j = 0; j < matr_res.m_ncols; ++j) {
			matr_res.m_pData[i * matr_res.m_ncols + j] = m_pData[j * m_ncols + i];
		}
	}
	return matr_res;
}

template <typename MYTYPE> TMatrix<MYTYPE> TMatrix<MYTYPE>::precond_matr2() {
	TMatrix<MYTYPE> matr_res(m_ncols, m_nrows);
	for (auto i = 0; i < matr_res.m_nrows - 3; ++i) {
		for (auto j = 0; j < matr_res.m_ncols - 3; ++j) {
			if (i == j)
				matr_res.m_pData[i * matr_res.m_ncols + j] = m_pData[i * m_ncols + j];
			else
				matr_res.m_pData[i * matr_res.m_ncols + j] = MYTYPE(0);
		}
	}
	for (auto i = 0; i < matr_res.m_nrows; ++i) {
		for (auto j = matr_res.m_ncols - 3; j < matr_res.m_ncols; ++j) {
			matr_res.m_pData[i * matr_res.m_ncols + j] = m_pData[i * m_ncols + j];
		}
	}
	for (auto i = matr_res.m_nrows - 3; i < matr_res.m_nrows; ++i) {
		for (auto j = 0; j < matr_res.m_ncols; ++j) {
			matr_res.m_pData[i * matr_res.m_ncols + j] = m_pData[i * m_ncols + j];
		}
	}
	return matr_res;
}

template <typename MYTYPE> int TMatrix<MYTYPE>::get_nrows() const{
	return m_nrows;
}

template <typename MYTYPE> int TMatrix<MYTYPE>::get_ncols() const{
	return m_ncols;
}


template <typename MYTYPE> TMatrix<MYTYPE>& TMatrix<MYTYPE>::operator += (const TMatrix<MYTYPE>& matr) {
	add(matr);
	return *this;
}

template <typename MYTYPE> TMatrix<MYTYPE>& TMatrix<MYTYPE>::operator -= (const TMatrix<MYTYPE>& matr) {
	subtract(matr);
	return *this;
}

template <typename MYTYPE> TMatrix<MYTYPE> operator + (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	return sum(matr1, matr2);
}

template <typename MYTYPE> TMatrix<MYTYPE> operator - (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	return diff(matr1, matr2);
}

template <typename MYTYPE> TMatrix<MYTYPE> operator * (const TMatrix<MYTYPE>& matr1, const TMatrix<MYTYPE>& matr2) {
	return matrix_product(matr1, matr2);
}

template <typename MYTYPE> std::ostream& operator << (std::ostream& out, const TMatrix<MYTYPE>& matr) {
	for (auto i = 0; i < matr.m_nrows; ++i) {
		for (auto j = 0; j < matr.m_ncols; ++j) {
			out << std::setw(12) << matr.m_pData[i * matr.m_ncols + j] << " ";
		}
		out << "\n";
	}
	out << "\n";

	return out;
}

template <typename MYTYPE> TMatrix<MYTYPE> TMatrix<MYTYPE>::operator - () {
	TMatrix<MYTYPE> matr_res(m_nrows, m_ncols, m_pData);
	for (auto i = 0; i < matr_res.m_elem_num; ++i) {
		matr_res.m_pData[i] = -matr_res.m_pData[i];
	}
	return matr_res;
}

template <typename MYTYPE> MYTYPE& TMatrix<MYTYPE>::operator () (const int row_ind, const int col_ind) {
	//if (row_ind < m_dim && col_ind < m_dim) {
		return m_pData[row_ind * m_ncols + col_ind];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return MYTYPE(0);
	//}
}

template <typename MYTYPE> const MYTYPE& TMatrix<MYTYPE>::operator () (const int row_ind, const int col_ind) const {
	//if (row_ind < m_dim && col_ind < m_dim) {
		return m_pData[row_ind * m_ncols + col_ind];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return MYTYPE(0);
	//}
}

template <typename MYTYPE> MYTYPE& TMatrix<MYTYPE>::operator [] (const int index) {
	//if (index < m_dim) {
	return m_pData[index];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return m_pData[0];
	//}
}

template <typename MYTYPE> const MYTYPE& TMatrix<MYTYPE>::operator [] (const int index) const {
	//if (index < m_dim) {
	return m_pData[index];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return MYTYPE(0);
	//}
}

template <typename MYTYPE> TMatrix<MYTYPE>& TMatrix<MYTYPE>::operator = (const TMatrix<MYTYPE>& matr) {
	if (this != &matr) {
		m_nrows = matr.m_nrows;
		m_ncols = matr.m_ncols;
		m_elem_num = matr.m_elem_num;
		delete[] m_pData;
		m_pData = new MYTYPE[m_elem_num];
		for (auto i = 0; i < m_elem_num; ++i) {
			m_pData[i] = matr.m_pData[i];
		}
	}
	return *this;
}

template <typename MYTYPE> TMatrix<MYTYPE>& TMatrix<MYTYPE>::operator = (TMatrix<MYTYPE>&& matr) {
	close();
	m_nrows = matr.m_nrows;
	m_ncols = matr.m_ncols;
	m_elem_num = matr.m_elem_num;
	m_pData = matr.m_pData;
	matr.m_nrows = 0;
	matr.m_ncols = 0;
	matr.m_elem_num = 0;
	matr.m_pData = nullptr;
	return *this;
}

template <typename MYTYPE> TVector<MYTYPE> operator * (const TMatrix<MYTYPE>& matr, const TVector<MYTYPE>& vec) {
	TVector<MYTYPE> vec_res(matr.m_nrows);
	if (matr.m_ncols == vec.get_dim()) {
		for (auto i = 0; i < matr.m_nrows; ++i) {
			for (auto j = 0; j < matr.m_ncols; ++j) {
				vec_res[i] += matr(i, j) * vec[j];
			}
		}
	}
	else {
		std::cout << "Unsuitable vector and matrix dimensions! The result is a zero vector!\n";
	}
	return vec_res;
}

template <typename MYTYPE> TMatrix<MYTYPE> operator * (const MYTYPE numb, const TMatrix<MYTYPE>& matr) {
	TMatrix<MYTYPE> matr_res(matr);
	for (auto i = 0; i < matr_res.m_elem_num; ++i) {
		matr_res.m_pData[i] *= numb;
	}
	return matr_res;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::set_as_column(const TVector<MYTYPE>& vec, int num_col) {
	if (m_nrows == vec.get_dim()) {
		if (num_col < m_ncols) {
			for (auto i = 0; i < m_nrows; ++i) {
				m_pData[i * m_ncols + num_col] = vec[i];
			}
		}
		else {
			std::cout << "The specified column number is greater than the number of matrix columns! It is not possible to set this vector as a matrix column!\n";
		}
	}
	else {
		std::cout << "It's impossible to set this vector as a column of the matrix: unsuitable vector and matrix dimensions! The original matrix is not changed!\n";
	}
}

template <typename MYTYPE> TVector<MYTYPE> TMatrix<MYTYPE>::get_column(int num_col) {
	TVector<MYTYPE> vec_res(m_nrows);
	for (auto i = 0; i < m_nrows; ++i) {
		vec_res[i] = m_pData[i * m_ncols + num_col];
	}
	return vec_res;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::remove_last_cols(int num_of_cols) {
	int new_ncols = m_ncols - num_of_cols;
	int new_elem_num = new_ncols * m_nrows;
	MYTYPE* new_pData = new MYTYPE[new_elem_num];
	for (auto j = 0; j < new_ncols; ++j) {
		for (auto i = 0; i < m_nrows; ++i) {
			new_pData[i * new_ncols + j] = m_pData[i * m_ncols + j];
		}
	}
	delete[] m_pData;
	m_pData = new_pData;
	new_pData = nullptr;
	m_ncols = new_ncols;
	m_elem_num = new_elem_num;
}

template <typename MYTYPE> void TMatrix<MYTYPE>::remove_last_rows(int num_of_rows) {
	int new_nrows = m_nrows - num_of_rows;
	int new_elem_num = m_ncols * new_nrows;
	MYTYPE* new_pData = new MYTYPE[new_elem_num];
	for (auto i = 0; i < new_nrows; ++i) {
		for (auto j = 0; j < m_ncols; ++j) {
			new_pData[i * m_ncols + j] = m_pData[i * m_ncols + j];
		}
	}
	delete[] m_pData;
	m_pData = new_pData;
	new_pData = nullptr;
	m_nrows = new_nrows;
	m_elem_num = new_elem_num;
}