#ifndef CPP8V_SQUARESPARSEMATRIX_H
#define CPP8V_SQUARESPARSEMATRIX_H
#include "SparseMatrix.h"
#include<cmath>
#include<exception>
#include<ostream>

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

    friend std::ostream& operator<<(std::ostream& os, SquareSparseMatrix & A) {
        std::map<size_t, std::map<size_t, double>> originalElementsA = A.getElements();
        std::map<size_t, std::map<size_t, double>> elementsA = {};

        size_t order = A.getOrder();

        elementsA = copyDoubleElements(A.getElements(), order);
        for (auto row : elementsA) {
            for (auto el : row.second) {
                os << el.second << "  ";
            }
            os << std::endl;
        }
        return os;

    }

    double calculateDeterminant() {
        return this->auxCalculateDeterminant(this->getElements());
    }

    SquareSparseMatrix<double> invertibleMatrix() {
        return SquareSparseMatrix<double>(auxInvertibleMatrix(this->getElements()));
    }


    friend SquareSparseMatrix<double> operator/(SquareSparseMatrix & A, SquareSparseMatrix & B) {
        if (A.getOrder() != B.getOrder()) {
            throw std::invalid_argument("Orders of A and B aren't equal");
        }
        std::map<size_t, std::map<size_t, double>> originalElementsA = A.getElements();
        std::map<size_t, std::map<size_t, double>> elementsA = {};

        size_t order = A.getOrder();

        elementsA = copyDoubleElements(A.getElements(), order);



        std::map<size_t, std::map<size_t, double>> elementsInvertibleB = \
        auxInvertibleMatrix(B.getElements());


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

    std::vector<double> solveEquation( std::vector<double> freeValues) {
        double determinant = calculateDeterminant();
        if (determinant == 0) {
            throw std::invalid_argument{"Irreversible matrix"};
        }
        size_t order = this->getOrder();
        std::map<size_t, std::map<size_t, double>> copiedElements = \
        copyDoubleElements(this->getElements(), order);

        size_t assistSize = order;
        for (size_t i = 0; i < order; i++) {
            copiedElements[i][assistSize] = freeValues[i];
        }
        // Here we make upper triangle matrix
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
                        break;
                    }
                    nonzeroIndex++;
                }
            }

            // Cycle for iterating over bottom rows, rows below i
            for (size_t j = i + 1; j < order; j++) {
                // Cycle for iterating over non-zero elements in a bottom row
                double multiplier = - copiedElements[j][i] / copiedElements[i][i];
                for (size_t k = i; k <= order; k++) {
                    copiedElements[j][k] += copiedElements[i][k] * multiplier;
                }
            }
        }
        // And here we make diagonal matrix
        for (size_t i = order - 1; i > 0; i--) {
            // Cycle for iterating over bottom rows, rows below i
            for (int j = (int)i - 1; j >= 0; j--) {
                // Cycle for iterating over non-zero elements in a bottom row
                double multiplier = - copiedElements[j][i] / copiedElements[i][i];
                copiedElements[j][i] += copiedElements[i][i] * multiplier;
                copiedElements[j][order] += copiedElements[i][order] * multiplier;
            }
        }
        //
        for (size_t i = 0; i < order; i++) {
            copiedElements[i][order] /= copiedElements[i][i];
        }
        //making vector of solves
        std::vector<T> Solves;
        for (size_t i = 0; i < this->getOrder(); i++) {
            Solves.push_back(copiedElements[i][order]);
        }
        return Solves;
    }



private:
    // This method is only used with dictionaries, not matrix objects!
    template<class Q>
    static std::map<size_t, std::map<size_t, double>> auxInvertibleMatrix\
    (const std::map<size_t, std::map<size_t, Q>> & matrixElements) {
        // Check if it worth it
        double determinant = auxCalculateDeterminant(matrixElements);
        if (determinant == 0) {
            throw std::invalid_argument{"Irreversible matrix"};
        }

        // Create a copy of the matrix, but with double instead of various Q type
        size_t rowsCount = matrixElements.size();
        size_t columnsCount = 0;
        for (auto row : matrixElements) {
            columnsCount = std::max(row.second.size(), columnsCount);
        }

        if (rowsCount != columnsCount) {
            throw std::invalid_argument("The matrix isn't square!");
        }
        size_t order = rowsCount;
        std::map<size_t, std::map<size_t, double>> copiedElements = copyDoubleElements(matrixElements, order);

        // calculating Matrix Of Algebraic Additions
        std::map<size_t, std::map<size_t, double>> matrixOfAlgebraicAdditions = \
        auxCreateMatrixOfAlgebraicAdditions(copiedElements, order);
        // transposition
        std::map<size_t, std::map<size_t, double>> TransposedMatrix = \
        transpose(matrixOfAlgebraicAdditions, order);
        // finally divide matrix by determinant and return the result
        divideMatrixByNumber(TransposedMatrix, determinant, order);
        return TransposedMatrix;
    }

    template<class Q>
    static double auxCalculateDeterminant(const std::map<size_t, std::map<size_t, Q>> & matrixElements) {
        // Create a copy of the matrix, but with double instead of T type
        size_t rowsCount = matrixElements.size();
        size_t columnsCount = 0;
        for (auto row : matrixElements) {
            columnsCount = std::max(row.second.size(), columnsCount);
        }

        if (rowsCount != columnsCount) {
            throw std::invalid_argument("The matrix isn't square!");
        }
        size_t order = rowsCount;
        std::map<size_t, std::map<size_t, double>> copiedElements = copyDoubleElements(matrixElements, order);

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

    // A function for copying elements with their type changed to double
    template<class Q>
    static std::map<size_t, std::map<size_t, double>> copyDoubleElements\
    (const std::map<size_t, std::map<size_t, Q>> & matrixElements, const size_t order) {
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

    // A method for making a matrix of alg. additions for a given matrix; You mustn't skip zero elements in
    // matrixElements variable !
    static std::map<size_t, std::map<size_t, double>> auxCreateMatrixOfAlgebraicAdditions\
    (const std::map<size_t, std::map<size_t, double>> & matrixElements, const size_t order) {
        std::map<size_t, std::map<size_t, double>> copiedElements = copyDoubleElements(matrixElements, order);
        std::map<size_t, std::map<size_t, double>> matrixOfAlgebraicAdditions = {};
        for (size_t excludedRow = 0; excludedRow < order; excludedRow++) {
            matrixOfAlgebraicAdditions[excludedRow] = {};
            for (size_t excludedColumn = 0; excludedColumn < order; excludedColumn++) {
                std::map<size_t, std::map<size_t, T>> minor = {};
                size_t i = 0, j;
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
        return matrixOfAlgebraicAdditions;
    }

    // You mustn't skip zero elements in
    // matrixElements variable !
    template<class Q>
    static std::map<size_t, std::map<size_t, Q>> transpose\
    (const std::map<size_t, std::map<size_t, Q>> & matrixElements, const size_t order) {
        std::map<size_t, std::map<size_t, Q>> TransposedMatrix = {};
        for (size_t rowIndex = 0; rowIndex < order; rowIndex++) {
            TransposedMatrix[rowIndex] = {};
        }
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                TransposedMatrix[i][j] = matrixElements.at(j).at(i);
            }
        }
        return TransposedMatrix;
    }

    // You mustn't skip zero elements in
    // matrixElements variable !
    static std::map<size_t, std::map<size_t, double>> divideMatrixByNumber\
    (std::map<size_t, std::map<size_t, double>> & matrixElements, const double divider, const size_t order) {
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                matrixElements.at(i).at(j) /= divider;
            }
        }
        return matrixElements;
    }
};


#endif //CPP8V_SQUARESPARSEMATRIX_H
