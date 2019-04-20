#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
// #include <omp.h>
#include <mpi.h>
#include <climits>

void MPI_update(std::vector<double> & x, std::vector <std::vector <int> > recieveRows, std::vector <std::vector <int> > sendRows,
                std::vector<int> globToLoc) {
    // to identify the sending process
    MPI_Request *recieveRequests = new MPI_Request [recieveRows.size()];
    MPI_Request *sendRequests = new MPI_Request [sendRows.size()];

    double **recieveMsg = new double * [recieveRows.size()]; 
    double **sendMsg = new double * [sendRows.size()];

    for (int i = 0; i < recieveRows.size(); i++) {
        recieveMsg[i] = new double [recieveRows[i].size()]; // array of values we want to recieve from process i
        sendMsg[i] = new double [sendRows[i].size()]; // array of values we want to send to process i

        for (int j = 0; j < sendRows[i].size(); j++) { // sendRows[i][j] - global row number we want to send
            sendMsg[i][j] = x[globToLoc[sendRows[i][j]]]; // in x we have only local numeration
        }

        MPI_Isend(sendMsg[i], sendRows[i].size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &(sendRequests[i]));
        MPI_Irecv(recieveMsg[i], recieveRows[i].size(), MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &(recieveRequests[i]));
    }

    MPI_Waitall(sendRows.size(), sendRequests, MPI_STATUS_IGNORE);
    MPI_Waitall(recieveRows.size(), recieveRequests, MPI_STATUS_IGNORE);

    for (int i = 0; i < recieveRows.size(); i++) {
        for (int j = 0; j < recieveRows[i].size(); j++) {
            x[globToLoc[recieveRows[i][j]]] = recieveMsg[i][j];
        }
    }

    for (int i = 0; i < recieveRows.size(); i++) {
        delete [] sendMsg[i];
        delete [] recieveMsg[i];
    }

    delete [] sendMsg;
    delete [] recieveMsg;

    delete [] sendRequests;
    delete [] recieveRequests;

    return;
}

class MPI_CompressedSparseRowMatrix {
private:
    int rowsN, colsN; // matrix size

    int nProcs; // total number of processes

    std::vector<int> dataStartRow; // data start of i row in columnNumbersRow and elementsRow; 
                                   // the last element is the number of non-zero elements in matrix
    std::vector<int> globColumnNumbersRow; // global column numbers of non-zero elements
    std::vector<int> locColumnNumbersRow; // local column numbers of non-zero elements according to globToLoc
    std::vector<double> elementsRow; // all non-zero elements in a row

    std::vector<int> rows; // local row number -> global row number
    std::vector<int> globToLoc; // global column number -> local column number

    std::vector<std::vector <int> > recieveRows; // rows numbers for recieve; procN -> rows numbers from procN process
    std::vector<std::vector <int> > sendRows; // rows numbers for send; procN -> rows numbers from procN process

    std::vector<int> halo; // halo-rows for the process

    std::vector<int> part; // which process owns the particular global row
public:
    MPI_CompressedSparseRowMatrix(int rowsN, int colsN, int nProcs): rowsN(rowsN), colsN(colsN), nProcs(nProcs) {
        std::vector<int> tmp;

        for (int i = 0; i < nProcs; i++) {
            recieveRows.push_back(tmp);
            sendRows.push_back(tmp);
        }
    }

    int rowsHaloSize() { return rowsN+halo.size(); }

    int rowsSize() { return rowsN; } // get local number of rows

    int haloSize() { return halo.size(); }

    int globColsSize() { return colsN; } // get global number of columns

    int getProcessesNumnber() { return nProcs; }

    void copyPart(std::vector<int> partCopy) { part = partCopy; }

    std::vector<int> getPart() { return part; }

    void fillPart(int gridX, int gridY, int gridZ, int partsX, int partsY, int partsZ) {
        for (int i = 0; i < gridX*gridY*gridZ; i++) {
            part.push_back(-1);
        }

        int nProcs = partsX * partsY * partsZ;

        for (int procN = 0; procN < nProcs; procN++) {
            int procXCord = procN % partsX, // process coordinates
                procYCord = (procN % (partsX*partsY)) / partsX,
                procZCord = procN / (partsX*partsY);

            int beginX = (gridX/partsX)*procXCord + std::min(procXCord, gridX%partsX),
                beginY = (gridY/partsY)*procYCord + std::min(procYCord, gridY%partsY),
                beginZ = (gridZ/partsZ)*procZCord + std::min(procZCord, gridZ%partsZ);

            int nProcX = (gridX/partsX) + (procXCord<gridX%partsX),
                nProcY = (gridY/partsY) + (procYCord<gridY%partsY),
                nProcZ = (gridZ/partsZ) + (procZCord<gridZ%partsZ);

            int endX = beginX + nProcX,
                endY = beginY + nProcY,
                endZ = beginZ + nProcZ;

            for (int z = beginZ; z < endZ; z++) {
                for (int y = beginY; y < endY; y++) {
                    for (int x = beginX; x < endX; x++) {
                        int globRow = z*(gridX*gridY) + y*gridX + x;
                        part[globRow] = procN;
                    }
                }
            }
        }
    }

    // prints part
    void printPart() {
        std::cout << "Vector Part - which process owns the particular global row:" << std::endl;
        for (int i = 0; i < part.size(); i++) {
            std::cout << "(" << part[i] << ": " << i << ") " << std::endl;
        }
    }

    void printHalo() {
        std::cout << "Halo-rows for the process:" << std::endl;
        for (int i = 0; i < halo.size(); i++) {
            std::cout << "(" << halo[i] << ": " << i+rows.size() << ") " << std::endl;
        }
    }

    std::vector <std::vector <int> > getRecieveRows() { return recieveRows; }

    void printRecieveRows() {
        std::cout << "Rows numbers for recieve:" << std::endl;
        for (int i = 0; i < recieveRows.size(); i++) {
            std::cout << i << ": ";
            for (int j = 0; j < recieveRows[i].size(); j++) {
                std::cout << recieveRows[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::vector <std::vector <int> > getSendRows() { return sendRows; }

    void printSendRows() {
        std::cout << "Rows numbers for send:" << std::endl;
        for (int i = 0; i < sendRows.size(); i++) {
            std::cout << i << ": ";
            for (int j = 0; j < sendRows[i].size(); j++) {
                std::cout << sendRows[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // prints dataStartRow
    void printDataStartRow() {
        std::cout << "IA" << std::endl;
        for (int i = 0; i < dataStartRow.size(); i++) {
            std::cout << dataStartRow[i] << " ";
        }
        std::cout << std::endl;
    }

    // prints globColumnNumbersRow
    void printGlobalColumnNumbersRow() {
        std::cout << "global JA" << std::endl;
        for (int i = 0; i < globColumnNumbersRow.size(); i++) {
            std::cout << globColumnNumbersRow[i] << " ";
        }
        std::cout << std::endl;
    }

    // prints locColumnNumbersRow
    void printLocalColumnNumbersRow() {
        std::cout << "local JA" << std::endl;
        for (int i = 0; i < locColumnNumbersRow.size(); i++) {
            std::cout << locColumnNumbersRow[i] << " ";
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

    // prints the part of the global matrix with global rows numbers and global columns numbers
    void printGlobalNumeration() {
        std::cout << "Number of rows: " << rowsN << std::endl;
        std::cout << "Number of cols: " << colsN << std::endl;

        std::cout << "The part of the global matrix:" << std::endl;

        std::cout.precision(4);
        std::cout << std::setw(8) << " ";
        for (int j = 0; j < colsN; j++) {
            std::cout << std::setw(7) << j << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < rowsN; i++) {
            if ((i != 0) and (rows[i] != rows[i-1]+1)) {
                std::cout << std::setw(7) << "..." << std::endl;
            }

            std::cout << std::setw(7) << rows[i] << " ";
            int dataIndex = dataStartRow[i];
            for (int j = 0; j < colsN; j++) {
                if ((dataIndex != dataStartRow[i+1]) and (j == globColumnNumbersRow[dataIndex])) {
                    std::cout << std::setw(7) << elementsRow[dataIndex] << " ";
                    dataIndex++;
                } else {
                    std::cout << std::setw(7) << 0 << " ";
                }
            }
            
            std::cout << std::endl;
        }
    }

    std::vector<int> getGlobToLoc() { return globToLoc; } 

    void pushBackRows(int globRow) {
        rows.push_back(globRow);
    }

    int getGlobRowNumber(int locRow) {
        return rows[locRow];
    }

    void fillCommunicationRows() { // fills recieveRows and sendRows
        // fills recieveRows and sendRows
        for (int locRow = 0; locRow < rowsN; locRow++) {
            for (int dataIndex = dataStartRow[locRow]; dataIndex < dataStartRow.at(locRow+1); dataIndex++) {
                int globCol = globColumnNumbersRow.at(dataIndex);
                 
                if (globToLoc[globCol] == -1) {
                    // this is halocell
                    recieveRows.at(part[globCol]).push_back(globCol);
                    sendRows.at(part[globCol]).push_back(rows[locRow]);
                }
            }
        }

        // deletes dublicates and sort cells (rows) to send and receive
        for (int i = 0; i < recieveRows.size(); i++) {
            std::set <int> tmpRecieve(recieveRows[i].begin(), recieveRows[i].end());
            recieveRows[i].assign(tmpRecieve.begin(), tmpRecieve.end());
            std::sort(recieveRows[i].begin(), recieveRows[i].end());

            std::set <int> tmpSend(sendRows[i].begin(), sendRows[i].end());
            sendRows[i].assign(tmpSend.begin(), tmpSend.end());  
            std::sort(sendRows[i].begin(), sendRows[i].end());         
        }
    }

    void fillHaloAndGlobToLoc() {
        // fills globToLoc with known rows
        for (int i = 0; i < colsN; i++) {
            globToLoc.push_back(-1);
        }
        for (int i = 0; i < rowsN; i++) {
            globToLoc[rows[i]] = i;
        }

        fillCommunicationRows();

        // forms halo and fills globToLoc to the end
        for (int i = 0; i < recieveRows.size(); i++) {
            for (int j = 0; j < recieveRows[i].size(); j++) {
                halo.push_back(recieveRows[i][j]);
                globToLoc[halo.back()] = rowsN + halo.size() - 1;
            }
        }
    }

    void globToLocColumnNumber() {
        fillHaloAndGlobToLoc();

        for (int i = 0; i < globColumnNumbersRow.size(); i++) {
            locColumnNumbersRow.push_back(globToLoc[globColumnNumbersRow[i]]);
        }
    }

    // changes the value of a non-zero element which exactly exists (in local rows and global columns)
    void changeValue(int locRow, int globCol, double element) {
        // std::cout << locRow << " " << globCol << " " << element << std::endl;
        for (int dataIndex = dataStartRow[locRow]; ; dataIndex++) {
            if (globCol == globColumnNumbersRow.at(dataIndex)) { // there is already an element at this place
                elementsRow[dataIndex] = element;
                return;
            }
        }
        std::cout << "There is no element with local row = " << locRow << ", global column = " << globCol << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // adds an element to the end during the sequential filling of the matrix (in local rows and global columns)
    double pushBack(int locRow, int globCol, double element) {
        // std::cout << locRow << " " << globCol << " " << element << std::endl;
        try { dataStartRow.at(locRow); }
        catch (...) { dataStartRow.push_back(elementsRow.size()); }

        globColumnNumbersRow.push_back(globCol);
        elementsRow.push_back(element);
        return element;
    }

    // (in local rows and global columns; element exactly exists!)
    double pushBackSinRowColumn(int locRow, int globRow, int globCol) {
        return pushBack(locRow, globCol, sin(globRow+globCol+1));
    }

    // (in local rows and global columns; element exactly exists!)
    double getElementInGlobalColumns(int locRow, int globCol) {
        for (int dataIndex = dataStartRow[locRow]; ; dataIndex++) {
            if (globCol == globColumnNumbersRow.at(dataIndex)) {
                return elementsRow[dataIndex];
            }
        }
        std::cout << "There is no element with local row = " << locRow << ", global column = " << globCol << ".\n";
        MPI_Finalize();
        exit(0);
    }    

    // (in local rows and local columns)
    double getElementInLocalColumns(int locRow, int locCol) {
        for (int dataIndex = dataStartRow[locRow]; ; dataIndex++) {
            if (locCol == locColumnNumbersRow.at(dataIndex)) {
                return elementsRow[dataIndex];
            }
        }
        std::cout << "There is no element with local row = " << locRow << ", local column = " << locCol << ".\n";
        MPI_Finalize();
        exit(0);
    }

    void addNumberOfElements() {
        dataStartRow.push_back(elementsRow.size());
    }

    double getNumberOfElements() {
        return elementsRow.size();
    }

    std::vector<double> & MPI_matrixVectorProduct(std::vector<double> & x, std::vector<double> & y) {
        // #pragma omp parallel for
        for (int i = 0; i < rowsSize(); i++) {
            y[i] = 0;
            for (int dataIndex = dataStartRow[i]; dataIndex < dataStartRow[i+1]; dataIndex++) {
                int j = locColumnNumbersRow[dataIndex];
                y[i] += elementsRow[dataIndex] * x[j];
            }
        }
        return y;
    }
};