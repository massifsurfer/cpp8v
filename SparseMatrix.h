#ifndef CPP8V_SPARSEMATRIX_H
#define CPP8V_SPARSEMATRIX_H

#include <iostream>
#include <map>
#include <cmath>
#include <random>
#include <exception>




template<class T>
class SparseMatrix
{
    size_t _rows;
    size_t _columns;
    std::map<size_t, std::map<size_t, T>> _elements;
public:
    SparseMatrix(size_t rows, size_t columns): _rows(rows), _columns(columns) {
        this->_elements = {};
        for (size_t i = 0; i < this->_rows; i++) {
            for (size_t j = 0; j < this->_columns; j++) {
                T f = (T)std::rand() % 11;
                if (f == 0)
                    continue;
                this->_elements[i][j] = f;
            }
        }
    }

    SparseMatrix(const std::map<size_t, std::map<size_t, T>> & elements) {
        size_t rowsCount = elements.size();
        size_t columnsCount = 0;
        for (auto row : elements) {
            columnsCount = std::max(row.size(), columnsCount);
        }
        this->_rows = rowsCount;
        this->_columns = columnsCount;
        this->_elements = elements;
    }

    size_t getRows() const {
        return this->_rows;
    }

    size_t getColumns() const {
        return this->_columns;
    }

    const std::map<size_t, std::map<size_t, T>> & getElements() {
        return this->_elements;
    }
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