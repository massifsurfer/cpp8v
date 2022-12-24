#ifndef CPP8V_SQUARESPARSEMATRIX_H
#define CPP8V_SQUARESPARSEMATRIX_H
#include "SparseMatrix.h"
#include<cmath>
#include<exception>

template<class T>
void printSparseMatrix(const std::map<size_t, std::map<size_t, T>> & elements) {

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
    SquareSparseMatrix(const std::map<size_t, std::map<size_t, T>> & elements): SparseMatrix<T>(elements) {
        if (this->getRows() != this->getColumns()) {
            throw std::invalid_argument("Only square matrix are acceptable");
        }
    }

    size_t getOrder() const {
        return this->getRows();
    }

    double calculateDeterminant() {
        return this->auxCalculateDeterminant(this->getElements());
    }

    SquareSparseMatrix<double> invertibleMatrix() {
        return SquareSparseMatrix<double>(auxInvertibleMatrix(this->getElements()));;
    }


    friend SquareSparseMatrix<double> operator/(SquareSparseMatrix & A, SquareSparseMatrix & B) {
        if (A.getOrder() != B.getOrder()) {
            throw std::invalid_argument("Orders of A and B aren't equal");
        }
        std::map<size_t, std::map<size_t, double>> originalElementsA = A.getElements();
        std::map<size_t, std::map<size_t, double>> elementsA = {};

        size_t order = A.getOrder();

        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            elementsA[rowIndex] = {};
            for (size_t columnIndex = 0; columnIndex < order; columnIndex++) {
                if (originalElementsA.at(rowIndex).find(columnIndex) == originalElementsA.at(rowIndex).end()) {
                    elementsA[rowIndex][columnIndex] = 0;
                }
                else {
                    elementsA[rowIndex][columnIndex] = (double) originalElementsA.at(rowIndex).at(columnIndex);
                }
            }
        }

        printSparseMatrix(elementsA);


        std::map<size_t, std::map<size_t, double>> elementsInvertibleB = auxInvertibleMatrix(B.getElements());

        printSparseMatrix(elementsInvertibleB);
        std::map<size_t, std::map<size_t, double>> resultElements = {};
        for (size_t rowIndex = 0; rowIndex < A.getOrder(); rowIndex++) {
            resultElements[rowIndex] = {};
        }
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                for (size_t k = 0; k < order; k++) {
                    resultElements[i][j] += elementsA[i][k]*elementsInvertibleB[k][j];
                }
            }
        }
        return resultElements;
    }


private:
    template<class K>
    static std::map<size_t, std::map<size_t, double>> getDoubleCopy(const std::map<size_t, std::map<size_t, K>> & matrixElements) {

        size_t rowsCount = matrixElements.size();
        size_t columnsCount = 0;
        for (auto row : matrixElements) {
            columnsCount = std::max(row.second.size(), columnsCount);
        }

        if (rowsCount != columnsCount) {
            throw std::invalid_argument("The matrix isn't square!");
        }
        size_t order = rowsCount;


        std::map<size_t, std::map<size_t, double>> copiedElements = {};
        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            copiedElements[rowIndex] = {};
            for (size_t columnIndex = 0; columnIndex < order; columnIndex++) {
                if (matrixElements.at(rowIndex).find(columnIndex) == matrixElements.at(rowIndex).end()) {
                    copiedElements[rowIndex][columnIndex] = 0;
                }
                else {
                    copiedElements[rowIndex][columnIndex] = (double) matrixElements.at(rowIndex).at(columnIndex);
                }
            }
        }
        return copiedElements;
    }



    // This method is only used with dictionaries, not matrix objects!
    static std::map<size_t, std::map<size_t, double>> auxInvertibleMatrix(const std::map<size_t, std::map<size_t, T>> & matrixElements) {

        // Check if it worth it
        double determinant = auxCalculateDeterminant(matrixElements);
        if (determinant == 0) {
            throw std::invalid_argument{"Irreversible matrix"};
        }

        // Create a copy of the matrix, but with double instead of T type
        std::map<size_t, std::map<size_t, double>> copiedElements = {};

        size_t rowsCount = matrixElements.size();
        size_t columnsCount = 0;
        for (auto row : matrixElements) {
            columnsCount = std::max(row.second.size(), columnsCount);
        }

        if (rowsCount != columnsCount) {
            throw std::invalid_argument("The matrix isn't square!");
        }
        size_t order = rowsCount;

        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            copiedElements[rowIndex] = {};
            for (size_t columnIndex = 0; columnIndex < order; columnIndex++) {
                if (matrixElements.at(rowIndex).find(columnIndex) == matrixElements.at(rowIndex).end()) {
                    copiedElements[rowIndex][columnIndex] = 0;
                }
                else {
                    copiedElements[rowIndex][columnIndex] = (double) matrixElements.at(rowIndex).at(columnIndex);
                }
            }
        }
        // calculating Matrix Of Algebraic Additions
        std::map<size_t, std::map<size_t, double>> matrixOfAlgebraicAdditions = {};
        for (size_t excludedRow = 0; excludedRow < order; excludedRow++) {
            matrixOfAlgebraicAdditions[excludedRow] = {};
            for (size_t excludedColumn = 0; excludedColumn < order; excludedColumn++) {
                std::map<size_t, std::map<size_t, T>> minor = {};
                size_t i = 0, j = 0;
                while (i < order) {
                    // Skip excluded row
                    if (i != excludedRow) {
                        if (i >  excludedRow) {
                            minor[i - 1] = {};
                        } else {
                            minor[i] = {};
                        }
                        j = 0;
                        while (j < order) {
                            // Skip excluded column and consider this fact in minor's indexing
                            if (i < excludedRow && j < excludedColumn) {
                                minor[i][j] = copiedElements[i][j];
                            }
                            else if (i > excludedRow && j < excludedColumn) {
                                minor[i - 1][j] = copiedElements[i][j];
                            }
                            else if(i < excludedRow && j > excludedColumn) {
                                minor[i][j - 1] = copiedElements[i][j];
                            }
                            else if(i > excludedRow && j > excludedColumn) {
                                minor[i - 1][j - 1] = copiedElements[i][j];
                            }
                            j++;
                        }
                    }
                    i++;
                }
                matrixOfAlgebraicAdditions[excludedRow][excludedColumn] = \
                auxCalculateDeterminant(minor) * std::pow(-1, excludedRow + excludedColumn);
            }
        }

        // transposition
        std::map<size_t, std::map<size_t, double>> TransposedMatrix = {};
        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            TransposedMatrix[rowIndex] = {};
        }
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                TransposedMatrix[i][j] = matrixOfAlgebraicAdditions[j][i] / determinant;
            }
        }
        return TransposedMatrix;
    }

    static double auxCalculateDeterminant(const std::map<size_t, std::map<size_t, T>> & matrixElements) {
        // Create a copy of the matrix, but with double instead of T type
        std::map<size_t, std::map<size_t, double>> copiedElements = {};

        size_t rowsCount = matrixElements.size();
        size_t columnsCount = 0;
        for (auto row : matrixElements) {
            columnsCount = std::max(row.second.size(), columnsCount);
        }

        if (rowsCount != columnsCount) {
            throw std::invalid_argument("The matrix isn't square!");
        }
        size_t order = rowsCount;

        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            copiedElements[rowIndex] = {};
            for (size_t columnIndex = 0; columnIndex < order; columnIndex++) {
                if (matrixElements.at(rowIndex).find(columnIndex) == matrixElements.at(rowIndex).end()) {
                    copiedElements[rowIndex][columnIndex] = 0;
                }
                else {
                    copiedElements[rowIndex][columnIndex] = (double) matrixElements.at(rowIndex).at(columnIndex);
                }
            }
        }


        double result = 1;
        // Iterate over all the rows for first step of Gaussian algorithm
        for (size_t i = 0; i < order - 1; i++) {

            // In this part we seek for a string with a non-zero first element. If there is a zero column, then the
            // determinant equals zero
            size_t nonzeroIndex = i;

            // If the pivot in the current row equals zero then start seeking for non-zero lower in the same column
            if (std::abs(copiedElements[nonzeroIndex][i]) < EPS) {

                nonzeroIndex++;
                while (nonzeroIndex < order) {
                    if (std::abs(copiedElements[nonzeroIndex][i]) > EPS) {
                        std::swap(copiedElements[i], copiedElements[nonzeroIndex]);
                        result *= -1;
                        break;
                    }
                    nonzeroIndex++;
                    if (nonzeroIndex == order)
                        return 0;
                }
            }

            // A cycle for iterating over bottom rows, rows below i
            for (size_t j = i + 1; j < order; j++) {
                // A cycle for iterating over non-zero elements in a bottom row
                double multiplier = - copiedElements[j][i] / copiedElements[i][i];
                for (size_t k = i; k < order; k++) {
                    copiedElements[j][k] += copiedElements[i][k] * multiplier;
                }
            }
        }

        // Finally multiply diagonal elements
        for (size_t i = 0; i < order; i++) {
            result *= copiedElements[i][i];
        }

        return result;
    }




};


#endif //CPP8V_SQUARESPARSEMATRIX_H
