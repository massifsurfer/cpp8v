#include <iostream>
#include "SquareSparseMatrix.h"

int main() {
    std::map<size_t, std::map<size_t, double>> A = {};
    A[0][0] = 1;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 2;
    SquareSparseMatrix<double> AS(A);

    std::map<size_t, std::map<size_t, double>> B = {};
    B[0][0] = 0;
    B[0][1] = 1;
    B[1][0] = 1;
    B[1][1] = 2;
    SquareSparseMatrix<double> BS(B);

    BS.solveEquation({1, 1});

    return 0;
}