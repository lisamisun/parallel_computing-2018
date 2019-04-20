#include <iostream>
#include <vector>
// #include <omp.h>

void MPI_testBasicOperationsTime(int nx, int ny, int nz, int px, int py, int pz, int procN) {
    MPI_RegularGridGenerator MPI_generatorA(nx, ny, nz, px, py, pz, procN);
    MPI_CompressedSparseRowMatrix * MPI_a = MPI_generatorA.MPI_generateDoubleCSRMatrix();

    int nTest = 10;
    double procT;
    double procTimeMPIDotProduct = 0, 
           procTimeMPILinearCombination = 0, 
           procTimeMPIMatrixVectorProduct = 0;
    double timeMPIDotProduct = 0, 
           timeMPILinearCombination = 0, 
           timeMPIMatrixVectorProduct = 0;

    // omp_set_num_threads(1);
    // test MPI parallel basic operations

    for (int i = 0; i < nTest; i++) {
        std::vector<double> x(MPI_a->rowsHaloSize(), 0), y(MPI_a->rowsHaloSize(), 0);
        for (int i = 0; i < MPI_a->rowsSize(); i++) {
            x[i] = sin(MPI_a->getGlobRowNumber(i));
            y[i] = cos(MPI_a->getGlobRowNumber(i));
        }

        procT = MPI_Wtime();
        for (int j = 0; j < nTest; j++) MPI_dotProduct(x, y, MPI_a->rowsSize(), MPI_a->haloSize());
        procTimeMPIDotProduct += MPI_Wtime() - procT;
        MPI_Reduce(&procTimeMPIDotProduct, &timeMPIDotProduct, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        procT = MPI_Wtime();
        for (int j = 0; j < nTest; j++) MPI_linearCombination(y, x, 1.00001, 0.99999, MPI_a->rowsSize(), MPI_a->haloSize());
        procTimeMPILinearCombination += MPI_Wtime() - procT;
        MPI_Reduce(&procTimeMPILinearCombination, &timeMPILinearCombination, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        procT = MPI_Wtime();
        for (int j = 0; j < nTest; j++) {
            MPI_update(x, MPI_a->getRecieveRows(), MPI_a->getSendRows(), MPI_a->getGlobToLoc()); // update x before testing
            MPI_a->MPI_matrixVectorProduct(x, y);
        }
        procTimeMPIMatrixVectorProduct += MPI_Wtime() - procT;
        MPI_Reduce(&procTimeMPIMatrixVectorProduct, &timeMPIMatrixVectorProduct, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (procN == 0) {
        std::cout << "testing MPI parallel ops:" << std::endl;
        std::cout.precision(3);
        std::cout << "axpy time=" << timeMPILinearCombination << "s" << std::endl;
        std::cout << "dot time=" << timeMPIDotProduct << "s" << std::endl;
        std::cout << "SpMV time=" << timeMPIMatrixVectorProduct << "s" << std::endl;
        std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}