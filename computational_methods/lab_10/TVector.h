#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
template <typename MYTYPE> class TVector 
{
protected:
    void init();
    void init(std::string filename);

public:
    int m_dim;
    MYTYPE* m_pData;
    void close();

	TVector() : m_dim{ 3 }, m_pData{ new MYTYPE[3] }{
		for (auto i = 0; i < m_dim; ++i) {
			m_pData[i] = MYTYPE(0);
		}
	}
	TVector(int dim, MYTYPE* pData);
	TVector(int dim);
	TVector(const TVector<MYTYPE>& vec);
	TVector(std::string filename);
	TVector(TVector<MYTYPE>&& vec);
	~TVector() {
		close();
	}

	void copy(int dim, MYTYPE* pData = nullptr);
	void copy(const TVector<MYTYPE>& vec);
	void output();
	double norm(char norm_type);
	void add(const TVector<MYTYPE>& vec);
	void subtract(const TVector<MYTYPE>& vec);
	void push_back(const MYTYPE elem);
	template <typename MYTYPE2> friend void swap(TVector<MYTYPE>& vec1, TVector<MYTYPE>& vec2);

	int get_dim() const;
	void zero_vector();
	void single_vector();

	void resize(int new_dim);

	template <typename MYTYPE2> friend MYTYPE scalar_product(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);
	template <typename MYTYPE2> friend TVector<MYTYPE> sum(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);
	template <typename MYTYPE2> friend TVector<MYTYPE> diff(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);

	TVector<MYTYPE>& operator += (const TVector<MYTYPE>& vec);
	TVector<MYTYPE>& operator -= (const TVector<MYTYPE>& vec);

	template <typename MYTYPE2> friend TVector<MYTYPE> operator + (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);
	template <typename MYTYPE2> friend TVector<MYTYPE> operator - (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);
	template <typename MYTYPE2> friend MYTYPE operator * (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2);
	template <typename MYTYPE2> friend TVector<MYTYPE> operator * (const MYTYPE numb, const TVector<MYTYPE>& vec);
	template <typename MYTYPE2> friend std::ostream& operator << (std::ostream& out, const TVector<MYTYPE>& vec);

	MYTYPE& operator [] (const int index);
	const MYTYPE& operator [] (const int index) const;

	TVector<MYTYPE> operator - ();
	TVector<MYTYPE>& operator = (const TVector<MYTYPE>& vec);
	TVector<MYTYPE>& operator = (TVector<MYTYPE>&& vec);

	void remove_last_elems(int num_of_elems);

	TVector(std::initializer_list<MYTYPE> lst) : m_pData(new MYTYPE[lst.size()]), m_dim(lst.size())
	{
		std::copy(lst.begin(), lst.end(), m_pData);
	}
};

template <typename MYTYPE> void TVector<MYTYPE>::init() {
	m_dim = 0;
	m_pData = nullptr;
}

template <typename MYTYPE> void TVector<MYTYPE>::init(std::string filename) {
	std::ifstream fin;
	fin.open(filename, std::ios::in);
	if (!fin.is_open())
	{
		std::cout << "Error file!\n";
	}
	else
	{
		fin >> m_dim;
		m_pData = new MYTYPE[m_dim];
		for (auto i = 0; i < m_dim; ++i) {
			fin >> m_pData[i];
		}
	}
	fin.close();
}

template <typename MYTYPE> void TVector<MYTYPE>::close() {
	if (m_dim) delete[] m_pData;
	init();
}

template <typename MYTYPE> void TVector<MYTYPE>::copy(int dim, MYTYPE* pData) {
	close();
	if (dim > 0) {
		m_dim = dim;
		m_pData = new MYTYPE[m_dim];
		if (pData) {
			for (auto i = 0; i < m_dim; ++i) {
				m_pData[i] = pData[i];
			}
		}
	}
}

template <typename MYTYPE> void TVector<MYTYPE>::copy(const TVector<MYTYPE>& vec) {
	copy(vec.m_dim, vec.m_pData);
}

template <typename MYTYPE> TVector<MYTYPE>::TVector(int dim, MYTYPE* pData) {
	init();
	copy(dim, pData);
}

template <typename MYTYPE> TVector<MYTYPE>::TVector(int dim) {
	m_dim = dim;
	m_pData = new MYTYPE[m_dim];
	for (auto i = 0; i < m_dim; ++i) {
		m_pData[i] = MYTYPE(0);
	}
}

template <typename MYTYPE> TVector<MYTYPE>::TVector(std::string filename) {
	init(filename);
}

template <typename MYTYPE> TVector<MYTYPE>::TVector(const TVector<MYTYPE>& vec) {
	init();
	copy(vec);
}

template <typename MYTYPE> TVector<MYTYPE>::TVector(TVector<MYTYPE>&& vec) {
	m_dim = vec.m_dim;
	m_pData = vec.m_pData;
	vec.m_dim = 0;
	vec.m_pData = nullptr;
}

template <typename MYTYPE> void TVector<MYTYPE>::output() {
	for (auto i = 0; i < m_dim; ++i) {
		std::cout << std::setw(12) << m_pData[i] << " ";
	}
	std::cout << "\n";
}

template <typename MYTYPE> MYTYPE scalar_product(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	MYTYPE scalar_product = MYTYPE(0);
	if (vec1.m_dim == vec2.m_dim) {
		for (auto i = 0; i < vec1.m_dim; ++i) {
			scalar_product += vec1.m_pData[i] * vec2.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of vectors! The scalar product := 0!";
	}
	return scalar_product;
}

template <typename MYTYPE> double TVector<MYTYPE>::norm(char norm_type) {
	MYTYPE norm = MYTYPE(0);
	switch (norm_type)
	{
	case '1':
		for (auto i = 0; i < m_dim; ++i) {
			norm += std::abs(m_pData[i]);
		}
		break;
	case '2':
		norm = sqrt(scalar_product(*this, *this));
		break;
	case '3':
		for (auto i = 0; i < m_dim; ++i) {
			if (std::abs(m_pData[i]) > norm) 
				norm = std::abs(m_pData[i]);
		}
		break;
	default:
		std::cout << "Wrong norm type! Norm := 0!";
		break;
	}
	return(norm);
}

template <typename MYTYPE> TVector<MYTYPE> sum(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	TVector<MYTYPE> vec_res(vec1);
	if (vec1.m_dim == vec2.m_dim) {
		for (auto i = 0; i < vec_res.m_dim; ++i) {
			vec_res.m_pData[i] += vec2.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of vectors! The result is an empty vector!\n";
		vec_res.close();
	}
	return vec_res;
}

template <typename MYTYPE> TVector<MYTYPE> diff(const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	TVector<MYTYPE> vec_res(vec1);
	if (vec1.m_dim == vec2.m_dim) {
		for (auto i = 0; i < vec_res.m_dim; ++i) {
			vec_res.m_pData[i] -= vec2.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of vectors! The result is an empty vector!\n";
		vec_res.close();
	}
	return vec_res;
}

template <typename MYTYPE> void TVector<MYTYPE>::add(const TVector<MYTYPE>& vec) {
	if (m_dim == vec.m_dim) {
		for (auto i = 0; i < m_dim; ++i) {
			m_pData[i] += vec.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of vectors! The original vector is not changed!";
	}
}

template <typename MYTYPE> void TVector<MYTYPE>::subtract(const TVector<MYTYPE>& vec) {
	if (m_dim == vec.m_dim) {
		for (auto i = 0; i < m_dim; ++i) {
			m_pData[i] -= vec.m_pData[i];
		}
	}
	else {
		std::cout << "Different dimensions of vectors! The original vector is not changed!";
	}
}

template <typename MYTYPE> void TVector<MYTYPE>::push_back(const MYTYPE elem) {
	int new_dim = m_dim + 1;
	MYTYPE* new_pData = new MYTYPE[new_dim];
	for (auto i = 0; i < m_dim; ++i) {
		new_pData[i] = m_pData[i];
	}
	new_pData[m_dim] = elem;
	delete[] m_pData;
	m_pData = new_pData;
	m_dim = new_dim;
	new_pData = nullptr;
}

template <typename MYTYPE> void swap(TVector<MYTYPE>& vec1, TVector<MYTYPE>& vec2) {
	std::swap(vec1.m_pData, vec2.m_pData);
}

template <typename MYTYPE> int TVector<MYTYPE>::get_dim() const{
	return m_dim;
}

template <typename MYTYPE> void TVector<MYTYPE>::zero_vector() {
	for (auto i = 0; i < m_dim; ++i) {
		m_pData[i] = MYTYPE(0);
	}
}

template <typename MYTYPE> void TVector<MYTYPE>::single_vector() {
	for (auto i = 0; i < m_dim; ++i) {
		m_pData[i] = MYTYPE(1);
	}
}

template <typename MYTYPE> void TVector<MYTYPE>::resize(int new_dim) {
	if (m_dim == new_dim) {
		return;
	}

	if (new_dim > 0) {
		delete[] m_pData;
	}

	m_dim = new_dim;
	if (m_dim > 0) {
		m_pData = new MYTYPE[m_dim];
		this->zero_vector();
	}
	else {
		m_pData = nullptr;
	}
}

template <typename MYTYPE> TVector<MYTYPE>& TVector<MYTYPE>::operator += (const TVector<MYTYPE>& vec) {
	add(vec);
	return *this;
}

template <typename MYTYPE> TVector<MYTYPE>& TVector<MYTYPE>::operator -= (const TVector<MYTYPE>& vec) {
	subtract(vec);
	return *this;
}

template <typename MYTYPE> TVector<MYTYPE> operator + (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	return sum(vec1, vec2);
}

template <typename MYTYPE> TVector<MYTYPE> operator - (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	return diff(vec1, vec2);
}

template <typename MYTYPE> MYTYPE operator * (const TVector<MYTYPE>& vec1, const TVector<MYTYPE>& vec2) {
	return scalar_product(vec1, vec2);
}

template <typename MYTYPE> TVector<MYTYPE> operator * (const MYTYPE numb, const TVector<MYTYPE>& vec) {
	TVector<MYTYPE> vec_res(vec);
	for (auto i = 0; i < vec_res.m_dim; ++i) {
		vec_res.m_pData[i] *= numb;
	}
	return vec_res;
}

template <typename MYTYPE> std::ostream& operator << (std::ostream& out, const TVector<MYTYPE>& vec) {
	for (auto i = 0; i < vec.m_dim; ++i) {
		out << std::setw(12) << vec.m_pData[i] << " ";
	}
	std::cout << "\n";
	return out;
}

template <typename MYTYPE> TVector<MYTYPE> TVector<MYTYPE>::operator - () {
	TVector<MYTYPE> vec_res(m_dim, m_pData);
	for (auto i = 0; i < vec_res.m_dim; ++i) {
		vec_res.m_pData[i] = -vec_res.m_pData[i];
	}
	return vec_res;
}

template <typename MYTYPE> MYTYPE& TVector<MYTYPE>::operator [] (const int index) {
	//if (index < m_dim) {
		return m_pData[index];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return m_pData[0];
	//}
}

template <typename MYTYPE> const MYTYPE& TVector<MYTYPE>::operator [] (const int index) const {
	//if (index < m_dim) {
		return m_pData[index];
	//}
	//else {
	//	std::cout << "Oversize!";
	//	return MYTYPE(0);
	//}
}

template <typename MYTYPE> TVector<MYTYPE>& TVector<MYTYPE>::operator = (const TVector<MYTYPE>& vec) {
	if (this != &vec) {
		m_dim = vec.m_dim;
		delete[] m_pData;
		m_pData = new MYTYPE[m_dim];
		for (auto i = 0; i < m_dim; ++i) {
			m_pData[i] = vec.m_pData[i];
		}
	}
	return *this;
}

template <typename MYTYPE> TVector<MYTYPE>& TVector<MYTYPE>::operator = (TVector<MYTYPE>&& vec) {
	close();
	m_dim = vec.m_dim;
	m_pData = vec.m_pData;
	vec.m_dim = 0;
	vec.m_pData = nullptr;
	return *this;
}

template <typename MYTYPE> void TVector<MYTYPE>::remove_last_elems(int num_of_elems) {
	int new_dim = m_dim - num_of_elems;
	MYTYPE* new_pData = new MYTYPE[new_dim];
	for (auto i = 0; i < new_dim; ++i) {
		new_pData[i] = m_pData[i];
	}
	delete[] m_pData;
	m_pData = new_pData;
	new_pData = nullptr;
	m_dim = new_dim;
}