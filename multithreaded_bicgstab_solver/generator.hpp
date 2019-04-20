#include "csr_matrix.hpp"
#include <cmath>

class DiagonalGenerator {
public:
    CompressedSparseRowMatrix * generateDiagonalMatrix(CompressedSparseRowMatrix * matrixA) {
        CompressedSparseRowMatrix * matrix = new CompressedSparseRowMatrix (matrixA->size());

        for (int i = 0; i < matrix->size(); i++) {
            matrix->pushBack(i, i, 1.0 / matrixA->getElement(i, i));
        }

        matrix->addNumberOfElements();
        return matrix;
    }
};

class RegularGridGenerator {
private:
    int gridX, gridY, gridZ;
public:
    RegularGridGenerator(int nx, int ny, int nz): gridX(nx), gridY(ny), gridZ(nz) {}

    CompressedSparseRowMatrix * generateDoubleCSRMatrix() {
        CompressedSparseRowMatrix * matrix = new CompressedSparseRowMatrix (gridX*gridY*gridZ);

        for (int z = 0; z < gridZ; z++) {
            for (int y = 0; y < gridY; y++) {
                for (int x = 0; x < gridX; x++) {
                    double diagonalElement = 0;

                    int row = z*(gridX*gridY) + y*gridX + x; // matrix row number
                    if (z > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row - gridX*gridY));
                    if (y > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row - gridX));
                    if (x > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row - 1));
                    matrix->pushBack(row, row, 0);
                    if (x < gridX-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row + 1));
                    if (y < gridY-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row + gridX));
                    if (z < gridZ-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, row + gridX*gridY));

                    diagonalElement *= 1.1;
                    matrix->changeValue(row, row, diagonalElement);
                }
            }
        }

        matrix->addNumberOfElements();
        return matrix;
    }
};