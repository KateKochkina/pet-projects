#pragma once

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>

template<typename T>
TVector<T>::TVector()
    : size(0), elem(nullptr) {
}

template<typename T>
TVector<T>::TVector(const std::string& file_name) {
    std::ifstream file(file_name);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    size_t n;
    if (!(file >> n)) {
        throw std::invalid_argument("Failed to read matrix size");
    }
    size = n;

    if (size == 0) {
        elem = nullptr;
    } else {
        elem = new T[size];
        for (size_t i = 0; i < size; ++i) {
            if (!(file >> elem[i])) {
                delete[] elem;
                throw std::invalid_argument("Failed to read vector elements");
            }
        }
    }

    file.close();
}

template<typename T>
TVector<T>::TVector(const TVector<T>& vector)
    : size(vector.size) {
    if (size == 0) {
        elem = nullptr;
    } else {
        elem = new T[size];
        for (size_t i = 0; i < size; ++i) {
            elem[i] = vector.elem[i];
        }
    }
}

template<typename T>
TVector<T>::TVector(std::initializer_list<T> lst) : elem(new T[lst.size()]), size(lst.size()) {
    if (size == 0) {
        elem = nullptr;
    } else {
        delete[] elem;
        elem = new T[size];
        for (size_t i = 0; i < size; ++i) {
            std::copy(lst.begin(), lst.end(), elem);
        }
    }
}

template<typename T>
TVector<T>::TVector(size_t size)
    : size(size) {
    if (size == 0) {
        elem = nullptr;
    } else {
        elem = new T[size];
    }
}

template<typename T>
TVector<T>::~TVector() {
    if (size > 0) {
        delete[] elem;
    }
}

template<typename T>
void TVector<T>::to_identity() {
    for (size_t i = 0; i < size; ++i) {
        elem[i] = 1;
    }
}

template<typename T>
size_t TVector<T>::get_size() const {
    return size;
}

template<typename T>
T TVector<T>::get_elem(size_t num) const {
    return elem[num];
}

template<typename T>
void TVector<T>::set_elem(size_t num, T new_elem) {
    elem[num] = new_elem;
}


template<typename T>
TVector<T>& TVector<T>::operator=(const TVector<T>& vector) {
    if (this != &vector) {
        if (vector.size != size) {
            delete[] elem;
        }

        size = vector.size;
        if (size == 0) {
            elem = nullptr;
        } else {
            delete[] elem;
            elem = new T[size];
            for (size_t i = 0; i < size; ++i) {
                elem[i] = vector.elem[i];
            }
        }
    }

    return *this;
}

template<typename T>
TVector<T>& TVector<T>::operator=(TVector<T>&& vector) noexcept {
    if (this != &vector) {
        if (size > 0) {
            delete[] elem;
        }

        elem = std::exchange(vector.elem, nullptr);
        size = std::exchange(vector.size, 0);
    }

    return *this;
}


template<typename T>
void TVector<T>::resize(size_t new_size) {
    if (size == new_size) {
        return;
    }

    if (size > 0) {
        delete[] elem;
    }

    size = new_size;
    if (size > 0) {
        elem = new T[size];
    } else {
        elem = nullptr;
    }
}

template<typename T>
void TVector<T>::to_zero() {
    for (size_t i = 0; i < size; ++i) {
        elem[i] = 0;
    }
}

template<typename T>
void TVector<T>::clear() {
    if (size > 0) {
        delete[] elem;
        elem = nullptr;
        size = 0;
    }
}

template<typename T>
T& TVector<T>::operator[](size_t index) {
    return elem[index];
}

template<typename T>
T& TVector<T>::operator[](size_t index) const {
    return elem[index];
}

template<typename T>
TVector<T> operator+(const TVector<T>& left, const TVector<T>& right) {
    if (left.size != right.size) {
        throw std::invalid_argument("size mismatch");
    }
    size_t size = left.size;

    TVector<T> sum_vector(size);
    for (size_t i = 0; i < size; ++i) {
        sum_vector.elem[i] = left.elem[i] + right.elem[i];
    }

    return sum_vector;
}

template<typename T>
TVector<T> operator-(const TVector<T>& left, const TVector<T>& right) {
    if (left.size != right.size) {
        throw std::invalid_argument("size mismatch");
    }
    size_t size = left.size;

    TVector<T> sum_vector(size);
    for (size_t i = 0; i < size; ++i) {
        sum_vector.elem[i] = left.elem[i] - right.elem[i];
    }

    return sum_vector;
}

//template <typename T>
//TVector<T> TVector<T>::operator*(const T num) {
//    TVector<T> mul_vector(size);
//    for (size_t i = 0; i < size; ++i) {
//        //mul_vector.elem[i] = elem[i] * num;
//    }
//    return mul_vector;
//}
template<typename T>
TVector<T> operator*(const TVector<T>& vector, const T num) {
    TVector<T> mul_vector(vector.get_size());
    for (size_t i = 0; i < vector.get_size(); ++i) {
        mul_vector.elem[i] = vector.elem[i] * num;
    }
    return mul_vector;
}

template<typename T>
TVector<T> TVector<T>::operator/(T num) {
    TVector<T> div_vector(size);

    for (size_t i = 0; i < size; ++i) {
        div_vector.elem[i] = elem[i] / num;
    }

    return div_vector;
}


template<typename T>
TVector<T> operator*(const T num, const TVector<T>& vector) {
    return vector * num;
}

template<typename T>
T operator*(const TVector<T>& left, const TVector<T>& right) {
    T sum = 0;
    for (size_t i = 0; i < left.get_size(); ++i) {
        sum += left[i] * right[i];
    }
    return sum;
}

template<typename T>
TVector<T> operator*(const TVector<TVector<T>>& matrix, const TVector<T>& vector) {
    size_t size = matrix.get_size();

    if (size != vector.get_size()) {
        throw std::invalid_argument("size mismatch");
    }

    TVector<T> mul_vector(size);
    for (size_t i = 0; i < size; ++i) {
        mul_vector[i] = 0;
        for (size_t k = 0; k < size; ++k) {
            mul_vector[i] += matrix[i][k] * vector[k];
        }
    }

    return mul_vector;
}

template<typename T>
T TVector<T>::norm(Norm_type norm_type) const {
    switch (norm_type) {
        case OCTAHEDRAL_NORM:
            return norm_octahedral();
        case SPHERICAL_NORM:
            return norm_spherical();
        case CUBIC_NORM:
            return norm_cubic();
        default:
            throw std::invalid_argument("Norm not found");
    }
}

template<typename T>
T TVector<T>::norm_octahedral() const {
    T norm = 0;
    for (size_t i = 0; i < size; ++i) {
        norm += std::abs(elem[i]);
    }
    return norm;
}

template<typename T>
T TVector<T>::norm_spherical() const {
    T norm = 0;
    for (size_t i = 0; i < size; ++i) {
        norm += elem[i] * elem[i];
    }
    norm = std::sqrt(norm);
    return norm;
}

template<typename T>
T TVector<T>::norm_cubic() const {
    T norm = 0;
    for (size_t i = 0; i < size; ++i) {
        if (norm < std::abs(elem[i])) {
            norm = std::abs(elem[i]);
        }
    }
    return norm;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const TVector<T>& vector) {
    for (size_t i = 0; i < vector.size; ++i) {
        out << std::setw(7) << vector.elem[i] << "\n";
    }

    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const TVector<TVector<T>>& matrix) {
    for (size_t i = 0; i < matrix.size; ++i) {
        for (size_t j = 0; j < matrix.size; ++j) {
            out << std::setw(7) << matrix.elem[i][j] << " ";
        }
        out << "\n";
    }

    return out;
}