#include <iostream>
#include "SquareSparseMatrix.h"

int main() {

    std::map<size_t, std::map<size_t, double>> A = {};
    A[0][0] = 1;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 1;

    std::map<size_t, std::map<size_t, double>> B = A;

    B = A;
    B[0][0] = -1;

    printSparseMatrix(B);

    printSparseMatrix(A);

    return 0;
}