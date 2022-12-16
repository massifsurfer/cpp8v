#include <iostream>
#include "SquareSparseMatrix.h"

int main() {

    SquareSparseMatrix<int> A(4);
    for (auto row : A.getElements()) {
        for (auto el : row.second) {
            std::cout << el.second << " ";
        }
        std::cout << std::endl;
    }
    std::cout << A.calculateDeterminant();
    return 0;
}