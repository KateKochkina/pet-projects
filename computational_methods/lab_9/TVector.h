#pragma once

#include<iostream>

enum Norm_type { OCTAHEDRAL_NORM, SPHERICAL_NORM, CUBIC_NORM };

template <typename T>
class TVector {
private:
    size_t size;
    T *elem;

public:
    //! init/release operations
    TVector();
    explicit TVector(const std::string &file_name);
    TVector(const TVector<T> &vector);
    TVector(std::initializer_list<T> lst);
    explicit TVector(size_t size);
    ~TVector();

    void resize(size_t new_size);
    void to_zero();
    void clear();
    void to_identity();

    //! basic operations
    size_t get_size() const;
    T get_elem(size_t num) const;
    void set_elem(size_t num, T new_elem);

    TVector<T> &operator=(const TVector<T> &vector);
    TVector<T> &operator=(TVector<T> &&vector) noexcept;

    T &operator[](size_t index);
    T &operator[](size_t index) const;

    template <typename U>
    friend TVector<U> operator+(const TVector<U> &left, const TVector<U> &right);
    
    template <typename U>
    friend TVector<U> operator-(const TVector<U> &left, const TVector<U> &right);

    template <typename U>
    friend TVector<U> operator*(const TVector<U> &vector, const U num);

    TVector<T> operator/(T num);

    template <typename U>
    friend TVector<U> operator*(const U num, const TVector<U> &vector);

    template <typename U>
    friend U operator*(const TVector<U> &left, const TVector<U> &right);

    template <typename U>
    friend TVector<U> operator*(const TVector<TVector<U>> &matrix, const TVector<U> &vector);

    T norm(Norm_type norm_type) const;
private:
    T norm_octahedral() const;
    T norm_spherical() const;
    T norm_cubic() const;

public:
    template <typename U>
    friend std::ostream &operator<<(std::ostream &out, const TVector<U> &vector);

    template <typename U>
    friend std::ostream &operator<<(std::ostream &out, const TVector<TVector<U>> &matrix);
};

#include "TVector.hpp"