#include "csr_matrix.hpp"
#include <cmath>

class MPI_DiagonalGenerator {
public:
    MPI_CompressedSparseRowMatrix * MPI_generateDiagonalMatrix(MPI_CompressedSparseRowMatrix * matrixA) {
        MPI_CompressedSparseRowMatrix * matrix = 
            new MPI_CompressedSparseRowMatrix (matrixA->rowsSize(), matrixA->globColsSize(), matrixA->getProcessesNumnber());

        matrix->copyPart(matrixA->getPart());

        for (int i = 0; i < matrix->rowsSize(); i++) { // fills the first ni rows (main part);
            int globRow = matrixA->getGlobRowNumber(i);
            matrix->pushBackRows(globRow);
            matrix->pushBack(i, globRow, 1.0 / matrixA->getElementInGlobalColumns(i, globRow));
        }

        matrix->addNumberOfElements();
        matrix->globToLocColumnNumber();

        return matrix;
    }
};

class MPI_RegularGridGenerator {
private:
    int gridX, gridY, gridZ; // regular grid dimension (? - нужно ли хранить)

    int procN, nProcs; // process number and total number of processes
    int partsX, partsY, partsZ; // number of parts for dividing the grid for each coordinate
    int procXCord, procYCord, procZCord; // process coordinates in processes grid 

    int nProcX, nProcY, nProcZ; // number of cells (matrix rows) in the process for each coordinate (? - нужно ли хранить)
public:
    MPI_RegularGridGenerator(int nx, int ny, int nz, int px, int py, int pz, int proc): 
        gridX(nx), gridY(ny), gridZ(nz), partsX(px), partsY(py), partsZ(pz), procN(proc) {
        nProcs = partsX * partsY * partsZ;

        procXCord = procN % partsX; // process coordinates
        procYCord = (procN % (partsX*partsY)) / partsX;
        procZCord = procN / (partsX*partsY);

        nProcX = (gridX/partsX) + (procXCord<gridX%partsX);
        nProcY = (gridY/partsY) + (procYCord<gridY%partsY);
        nProcZ = (gridZ/partsZ) + (procZCord<gridZ%partsZ);
    }

    MPI_CompressedSparseRowMatrix * MPI_generateDoubleCSRMatrix() {
        MPI_CompressedSparseRowMatrix * matrix = 
            new MPI_CompressedSparseRowMatrix(nProcX*nProcY*nProcZ, gridX*gridY*gridZ, nProcs);

        matrix->fillPart(gridX, gridY, gridZ, partsX, partsY, partsZ);

        int beginX = (gridX/partsX)*procXCord + std::min(procXCord, gridX%partsX),
            beginY = (gridY/partsY)*procYCord + std::min(procYCord, gridY%partsY),
            beginZ = (gridZ/partsZ)*procZCord + std::min(procZCord, gridZ%partsZ);

        int endX = beginX + nProcX,
            endY = beginY + nProcY,
            endZ = beginZ + nProcZ;

        int row = 0; // local row number

        for (int z = beginZ; z < endZ; z++) {
            for (int y = beginY; y < endY; y++) {
                for (int x = beginX; x < endX; x++) {
                    double diagonalElement = 0;

                    int globRow = z*(gridX*gridY) + y*gridX + x; // matrix row global number

                    matrix->pushBackRows(globRow);

                    if (z > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow - gridX*gridY));
                    if (y > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow - gridX));
                    if (x > 0) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow - 1));

                    matrix->pushBack(row, globRow, 0);

                    if (x < gridX-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow + 1));
                    if (y < gridY-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow + gridX));
                    if (z < gridZ-1) diagonalElement += fabs(matrix->pushBackSinRowColumn(row, globRow, globRow + gridX*gridY));

                    diagonalElement *= 1.1;

                    /* matrix->printDataStartRow();
                    matrix->printGlobalColumnNumbersRow();
                    matrix->printElementsRow();
                    */

                    matrix->changeValue(row, globRow, diagonalElement);

                    row++;
                }
            }
        }
        matrix->addNumberOfElements();
        matrix->globToLocColumnNumber();
        return matrix;
    }

};