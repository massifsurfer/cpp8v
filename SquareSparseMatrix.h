#ifndef CPP8V_SQUARESPARSEMATRIX_H
#define CPP8V_SQUARESPARSEMATRIX_H
#include "SparseMatrix.h"
#include<cmath>
#include<exception>

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

    std::map<size_t, std::map<size_t, double>> InvertibleMatrix() {

        double determinant = calculateDeterminant();
        if (determinant = 0) {
            throw std::invalid_argument{"Irreversible matrix"};
        }

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

        // searching Matrix Of Algebraic Additions
        std::map<size_t, std::map<size_t, double>> MatrixOfAlgebraicAdditions = {};
        std::map<size_t, std::map<size_t, double>> Minor = {};

        for (size_t i = 0; i < this->getOrder(); i++) {
            for (size_t j = 0; j < this->getOrder(); j++) {
                Minor = copiedElements;
                for (size_t k = 0; k < this->getOrder(); k++) {
                    Minor[k].erase(j);
                }
                Minor.erase(i);
                //MatrixOfAlgebraicAdditions[i][j] = std::pow(-1,i+j)*(Minor.calculateDeterminant());
            }
        }
        //transposition
        std::map<size_t, std::map<size_t, double>> TransposedMatrix = {};
        for (size_t rowIndex = 0; rowIndex < this->getOrder(); rowIndex++) {
            TransposedMatrix[rowIndex] = {};
        }
        for (size_t i = 0; i < this->getOrder(); i++) {
            for (size_t j = 0; j < this->getOrder(); j++) {
                TransposedMatrix[i][j] = MatrixOfAlgebraicAdditions[j][i];
            }
        }


        //last step
        for (size_t i = 0; i < this->getOrder(); i++) {
            for (size_t j = 0; j < this->getOrder(); j++) {
                MatrixOfAlgebraicAdditions[i][j] *= 1 / determinant;
            }
        }

        return MatrixOfAlgebraicAdditions;





    }

    friend SquareSparseMatrix operator*(SquareSparseMatrix A, SquareSparseMatrix B ) {
        if (A->getOrder() != B->getOrder()) {
            throw std::invalid_argument("Orders of A and B aren't equal");
        }
        B = B.InvertibleMatrix();
        SquareSparseMatrix C = {};
        for (size_t rowIndex = 0; rowIndex < A->getOrder(); rowIndex++) {
            C[rowIndex] = {};
        }
        for (size_t i = 0; i < A->getOrder(); i++) {
            for (size_t j = 0; j < A->getOrder(); j++) {
                for (size_t k = 0; k < A->getOrder(); k++) {
                    C[i][j] = A[i][k]*B[k][j];
                }
            }
        }
        return C;

    }
};


#endif //CPP8V_SQUARESPARSEMATRIX_H
