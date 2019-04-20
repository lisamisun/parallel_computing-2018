#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <omp.h>

class CompressedSparseRowMatrix {
private:
    int matrixSize; // matrix size
    std::vector<int> dataStartRow; // data start of i row in columnNumbersRow and elementsRow; the last element is the number of non-zero elements in matrix
    std::vector<int> columnNumbersRow; // column numbers of non-zero elements
    std::vector<double> elementsRow; // all non-zero elements in a row

public:
    CompressedSparseRowMatrix(int n): matrixSize(n) {}

    int size() {return matrixSize;}

    // prints dataStartRow
    void printDataStartRow() {
        std::cout << "IA" << std::endl;
        for (int i = 0; i < dataStartRow.size(); i++) {
            std::cout << dataStartRow[i] << " ";
        }
        std::cout << std::endl;
    }

    // prints columnNumbersRow
    void printColumnNumbersRow() {
        std::cout << "JA" << std::endl;
        for (int i = 0; i < columnNumbersRow.size(); i++) {
            std::cout << columnNumbersRow[i] << " ";
        }
        std::cout << std::endl;
    }

    // prints elementsRow
    void printElementsRow() {
        std::cout << "A" << std::endl;
        std::cout.precision(4);
        for (int i = 0; i < elementsRow.size(); i++) {
            std::cout << elementsRow[i] << " ";
        }
        std::cout << std::endl;
    }

    // prints the full matrix
    void print() {
        std::cout << "Matrix size: " << matrixSize << std::endl;
        std::cout << "Matrix:" << std::endl;

        std::cout.precision(4);
        std::cout << std::setw(8) << " ";
        for (int j = 0; j < matrixSize; j++) {
            std::cout << std::setw(7) << j << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < matrixSize; i++) {
            std::cout << std::setw(7) << i << " ";
            int dataIndex = dataStartRow[i];
            for (int j = 0; j < matrixSize; j++) {
                if ((dataIndex != dataStartRow[i+1]) and (j == columnNumbersRow[dataIndex])) {
                    std::cout << std::setw(7) << elementsRow[dataIndex] << " ";
                    dataIndex++;
                } else {
                    std::cout << std::setw(7) << 0 << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    // changes the value of a non-zero element which exactly exists
    void changeValue (int row, int column, double element) {
        for (int i = dataStartRow[row]; ; i++) {
            if (column == columnNumbersRow.at(i)) { // there is already an element at this place
                elementsRow[i] = element;
                return;
            }
        }
    }

    // adds an element to the end during the sequential filling of the matrix
    double pushBack(int row, int column, double element) {
        try { dataStartRow.at(row); }
        catch (...) { dataStartRow.push_back(elementsRow.size()); }

        columnNumbersRow.push_back(column);
        elementsRow.push_back(element);
        return element;
    }

    double pushBackSinRowColumn(int row, int column) {
        double a = pushBack(row, column, sin(row+column+1));
        return a;
    }

    double getElement(int row, int column) {
        for (int dataIndex = dataStartRow[row]; dataIndex < dataStartRow[row+1]; dataIndex++) {
            if (column == columnNumbersRow[dataIndex]) {
                return elementsRow[dataIndex];
            }
        }
        throw "There are no such elements.";
    }

    void addNumberOfElements () {
        dataStartRow.push_back(elementsRow.size());
    }

    double getNumberOfElements () {
        return elementsRow.size();
    }

    std::vector<double> & matrixVectorProduct(std::vector<double> & x, std::vector<double> & y) {
        for (int i = 0; i < matrixSize; i++) {
            y[i] = 0;
            for (int dataIndex = dataStartRow[i]; dataIndex < dataStartRow[i+1]; dataIndex++) {
                int j = columnNumbersRow[dataIndex];
                y[i] += elementsRow[dataIndex] * x[j];
            }
        }
        return y;
    }

    std::vector<double> & matrixVectorProduct_par(std::vector<double> & x, std::vector<double> & y) {
        #pragma omp parallel for
        for (int i = 0; i < matrixSize; i++) {
            y[i] = 0;
            for (int dataIndex = dataStartRow[i]; dataIndex < dataStartRow[i+1]; dataIndex++) {
                int j = columnNumbersRow[dataIndex];
                y[i] += elementsRow[dataIndex] * x[j];
            }
        }
        return y;
    }
};