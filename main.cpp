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
    //std::cout << A.calculateDeterminant();


    std::map<size_t, std::map<size_t, double>> test1 = {};
    for (size_t i = 0; i < 3; i++) {
        test1[i] = {};
        std::cout << '\n';
        for (size_t j = 0; j < 3; j++) {
            test1[i][j] = j;
            std::cout << test1[i][j];
        }
    }
    std::map<size_t, std::map<size_t, double>> TransposedMatrix = {};
    for (size_t rowIndex = 0; rowIndex < 3; rowIndex++) {
        TransposedMatrix[rowIndex] = {};
    }
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            TransposedMatrix[i][j] = test1[j][i];
        }
    }

    for (size_t i = 0; i < 3; i++) {
        std::cout << '\n';
        for (size_t j = 0; j < 3; j++) {
            std::cout << TransposedMatrix[i][j];
        }
    }

    return 0;
}