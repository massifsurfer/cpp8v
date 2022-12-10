#ifndef CPP8V_SPARSEMATRIX_H
#define CPP8V_SPARSEMATRIX_H

#include <iostream>
#include <map>


template<class T>
class SparseMatrix
{
    size_t _rows;
    size_t _columns;
    std::map<std::pair<size_t, size_t>, T> _elements;
public:
    SparseMatrix(size_t rows, size_t columns);
    // методы для доступа к матрице
    const T& at(size_t row, size_t column) const;
    T& at(size_t row, size_t column);
    SparseMatrix<T>& operator+=(const SparseMatrix<T>&);
    SparseMatrix<T>& operator-=(const SparseMatrix<T>&);
    SparseMatrix<T>& operator*=(const SparseMatrix<T>&);
    // проверка на равенство
    bool operator==(const SparseMatrix<T>&) const;
    bool operator!=(const SparseMatrix<T>&) const;
    SparseMatrix<T>& operator*=(const T&);
    SparseMatrix<T>& operator/=(const T&);
    // вывод в поток
    friend std::ostream& operator<<(std::ostream&, const SparseMatrix<T>&);
    // ранг матрицы
    size_t rang() const;
};

template<class T>
SparseMatrix<T>& operator+(const SparseMatrix<T>& left, const SparseMatrix<T>& right) {
    // TODO:
}

template<class T>
SparseMatrix<T>& operator-(const SparseMatrix<T>& left, const SparseMatrix<T>& right) {
    // TODO:
}

template<class T>
SparseMatrix<T>& operator*(const T& left, const SparseMatrix<T>& right) {
    // TODO:
}



#endif //CPP8V_SPARSEMATRIX_H
