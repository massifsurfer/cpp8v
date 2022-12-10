#ifndef CPP8V_SQUARESPARSEMATRIX_H
#define CPP8V_SQUARESPARSEMATRIX_H
#include "SparseMatrix.h"

template<class T>
class SquareSparseMatrix: public SparseMatrix<T>{
public:
    SquareSparseMatrix(size_t n): SparseMatrix<T>(n, n){};
};


#endif //CPP8V_SQUARESPARSEMATRIX_H
