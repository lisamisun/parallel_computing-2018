#include <iostream>
#include <vector>
#include <omp.h>

void testBasicOperationsSpeedup(int nx, int ny, int nz) {
    RegularGridGenerator generatorA(nx, ny, nz);
    CompressedSparseRowMatrix * a = generatorA.generateDoubleCSRMatrix();

    int nTest = 20, N = nx*ny*nz;
    double t;
    double tDotProductSeq = 0, 
           tLinearCombinationSeq = 0, 
           tMatrixVectorProductSeq = 0;

    const double dotProducFlop = nTest * nTest * N * 2 * 1E-9;
    const double linearCombinationFlop = nTest * nTest * N * 3 * 1E-9;
    const double matrixVectorProductFlop = nTest * nTest * a->getNumberOfElements() * 2 * 1E-9;

    // test sequential basic operations
    omp_set_num_threads(1);

    for (int i = 0; i < nTest; i++) {
        std::vector<double> x, y;
        for (int i = 0; i < nx*ny*nz; i++) {
            x.push_back(sin(i));
            y.push_back(cos(i));
        }

        t = omp_get_wtime();
        for (int j = 0; j < nTest; j++) dotProduct(x, y);
        tDotProductSeq += omp_get_wtime() - t;

        t = omp_get_wtime();
        for (int j = 0; j < nTest; j++) linearCombination(y, x, 1.00001, 0.99999);
        tLinearCombinationSeq += omp_get_wtime() - t;

        t = omp_get_wtime();
        for (int j = 0; j < nTest; j++) a->matrixVectorProduct(x, y);
        tMatrixVectorProductSeq += omp_get_wtime() - t;
    }

    std::cout << "testing sequential ops:" << std::endl;
    std::cout.precision(3);
    std::cout << "axpy time=" << tLinearCombinationSeq << "s GFLOPS=" << linearCombinationFlop/tLinearCombinationSeq << std::endl;
    std::cout << "dot time=" << tDotProductSeq << "s GFLOPS=" << dotProducFlop/tDotProductSeq << std::endl;
    std::cout << "SpMV time=" << tMatrixVectorProductSeq << "s GFLOPS=" << matrixVectorProductFlop/tMatrixVectorProductSeq << std::endl;
    std::cout << std::endl;


    // test parallel basic operations
    const int NTR = omp_get_num_procs();
    
    for (int ntr = 2; ntr <= NTR; ntr += 2) {
        std::vector<double> x, y;
        for (int i = 0; i < nx*ny*nz; i++) {
            x.push_back(sin(i));
            y.push_back(cos(i));
        }

        omp_set_dynamic(0);
        omp_set_num_threads(ntr);
        double tDotProductPar = 0, 
               tLinearCombinationPar = 0, 
               tMatrixVectorProductPar = 0;

        for (int i = 0; i < nTest; i++){
            t = omp_get_wtime();
            for (int j = 0; j < nTest; j++) linearCombination_par(y, x, 1.00001, 0.99999);
            tLinearCombinationPar += omp_get_wtime() - t;

            t = omp_get_wtime();
            for (int j = 0; j < nTest; j++) dotProduct_par(x, y);
            tDotProductPar += omp_get_wtime() - t;

            t = omp_get_wtime();
            for (int j = 0; j < nTest; j++) a->matrixVectorProduct_par(x, y);
            tMatrixVectorProductPar += omp_get_wtime() - t;
        }

        std::cout << "testing parallel ops for ntr=" << ntr << ":" << std::endl;
        std::cout.precision(3);
        std::cout << "axpy time=" << tLinearCombinationPar << "s GFLOPS=" << linearCombinationFlop/tLinearCombinationPar << " Speedup=" << tLinearCombinationSeq/tLinearCombinationPar << "X" << std::endl;
        std::cout << "dot time=" << tDotProductPar << "s GFLOPS=" << dotProducFlop/tDotProductPar << " Speedup=" << tDotProductSeq/tDotProductPar << "X" << std::endl;
        std::cout << "SpMV time=" << tMatrixVectorProductPar << "s GFLOPS=" << matrixVectorProductFlop/tMatrixVectorProductPar << " Speedup=" << tMatrixVectorProductSeq/tMatrixVectorProductPar << "X" << std::endl;
        std::cout << std::endl;
    }

}