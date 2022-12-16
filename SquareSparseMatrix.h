#ifndef CPP8V_SQUARESPARSEMATRIX_H
#define CPP8V_SQUARESPARSEMATRIX_H
#include "SparseMatrix.h"

template<class T>
void printSparseMatrix(std::map<size_t, std::map<size_t, T>> elements) {

    for (auto row : elements) {
        for (auto el : row.second) {
            std::cout << el.second << " ";
        }
        std::cout << std::endl;
    }
}

double EPS = 1.e-10;

template<class T>
class SquareSparseMatrix: public SparseMatrix<T>{
public:
    SquareSparseMatrix(size_t n): SparseMatrix<T>(n, n){};

    size_t getOrder() const {
        return this->getRows();
    }

    double calculateDeterminant() {
        // Create a copy of the matrix, but with double instead of T type
        std::map<size_t, std::map<size_t, double>> copiedElements = {};

        auto elements = this->getElements();
        for (size_t rowIndex = 0; rowIndex < this->getOrder(); rowIndex++) {
            copiedElements[rowIndex] = {};
            for (size_t columnIndex = 0; columnIndex < this->getOrder(); columnIndex++) {
                if (elements[rowIndex].find(columnIndex) == elements[rowIndex].end()) {
                    copiedElements[rowIndex][columnIndex] = 0;
                }
                else {
                    copiedElements[rowIndex][columnIndex] = (double) elements[rowIndex][columnIndex];
                }
            }
        }


        double result = 1;
        // Iterate over all the rows for first step of Gaussian algorithm
        for (size_t i = 0; i < this->getOrder() - 1; i++) {

            // In this part we seek for a string with a non-zero first element. If there is a zero column, then the
            // determinant equals zero
            size_t nonzeroIndex = i;

            // If the pivot in the current row equals zero then start seeking for non-zero lower in the same column
            if (std::abs(copiedElements[nonzeroIndex][i]) < EPS) {

                nonzeroIndex++;
                while (nonzeroIndex < this->getOrder()) {
                    if (std::abs(copiedElements[nonzeroIndex][i]) > EPS) {
                        std::swap(copiedElements[i], copiedElements[nonzeroIndex]);
                        result *= -1;
                        break;
                    }
                    nonzeroIndex++;
                    if (nonzeroIndex == this->getOrder())
                        return 0;
                }
            }

            // Cycle for iterating over bottom rows, rows below i
            for (size_t j = i + 1; j < this->getOrder(); j++) {
                // Cycle for iterating over non-zero elements in a bottom row
                double multiplier = - copiedElements[j][i] / copiedElements[i][i];
                for (size_t k = i; k < this->getOrder(); k++) {
                    copiedElements[j][k] += copiedElements[i][k] * multiplier;
                }
            }
        }

        // Finally multiply diagonal elements
        for (size_t i = 0; i < this->getOrder(); i++) {
            result *= copiedElements[i][i];
        }

        return result;
    }
};


#endif //CPP8V_SQUARESPARSEMATRIX_H
