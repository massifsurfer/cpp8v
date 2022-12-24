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
    std::vector<double> S = { 1, 3};
    std::vector<double> S1 = C.SLAU(S);
    for (int i = 0; i < 2; i++) {
        std::cout<<S1[i]<< '\n';
    }
    return 0;
}