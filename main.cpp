#include <iostream>
#include "SquareSparseMatrix.h"

int main() {
    std::map<size_t, std::map<size_t, double>> B = {};
    B[0][0] = 1;
    B[0][1] = 1;
    B[1][0] = 1;
    B[1][1] = 2;
    SquareSparseMatrix<double> C(B);
    printSparseMatrix(C.invertibleMatrix().getElements());
    return 0;
}