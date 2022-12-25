#include <iostream>
#include "SquareSparseMatrix.h"
#include<vector>
#include<cassert>


int main() {
    std::map<size_t, std::map<size_t, double>> A = {};
    A[0][0] = 1;
    A[0][1] = 1;
    A[0][2] = 0;
    A[1][0] = 1;
    A[1][1] = 2;
    A[1][2] = 1;
    A[2][0] = 0;
    A[2][1] = 0;
    A[2][2] = 1;
    SquareSparseMatrix<double> AS(A);

    std::map<size_t, std::map<size_t, double>> B = {};
    B[0][0] = 2;
    B[0][1] = 0;
    B[0][2] = 1;
    B[1][0] = 0;
    B[1][1] = 2;
    B[1][2] = 0;
    B[2][0] = 1;
    B[2][1] = 1;
    B[2][2] = 1;
    SquareSparseMatrix<double> BS(B);

    std::vector<double> SolvedAS = AS.solveEquation({1, 2, 1});
    std::cout << "Matrix A: " << std::endl;
    std::cout << AS;
    std::cout << "Matrix B: " << std::endl;
    std::cout << BS;
    std::cout << '\n' << "Determinant of A: " << std::endl;
    std::cout << AS.calculateDeterminant();
    SquareSparseMatrix<double> invAS = AS.invertibleMatrix();
    std::cout << '\n' << "Inverse Matrix A: " << std::endl;
    std::cout << invAS;
    std::cout << '\n' << "Solves of A: " << std::endl;
    for(size_t i = 0; i < 3; i++) {
        std::cout << SolvedAS[i] << " ";
    }
    std::cout << '\n' << "The result of dividing matrix A by matrix B: " << std::endl;
    SquareSparseMatrix<double> CS = AS/BS;
    std::cout << CS;
    return 0;
}