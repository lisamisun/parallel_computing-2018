#include <stdlib.h>
#include <stdio.h>

#include "parse_args.hpp"
#include "lin_algebra.hpp"
#include "generator.hpp"
#include "test_basic_operations.hpp"

int solverBiCGSTAB(std::vector<double> & solutionVector, int sizeN, CompressedSparseRowMatrix * matrixA, std::vector<double> & vectorBB, double criterionTol, int iterationsMax) {
    DiagonalGenerator generatorDD;
    CompressedSparseRowMatrix * dd = generatorDD.generateDiagonalMatrix(matrixA);
    
    std::vector<double> pp(sizeN, 0), pp2(sizeN, 0), rr = vectorBB, rr2 = vectorBB, tt(sizeN, 0), vv(sizeN, 0), ss(sizeN, 0), ss2(sizeN, 0);

    double initres = sqrt(dotProduct(vectorBB, vectorBB)), res = initres;
    double mineps = 1E-15;
    double eps = std::max(criterionTol*initres, mineps);

    double Rhoi_1 = 1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, Rhoi_2 = 1.0, alphai_1 = 1.0, wi_1 = 1.0;
    double RhoMin = 1E-60;

    int info = 1;

    std::cout << "Solver_BiCGSTAB: " << ": initres: " << initres << "; eps: " << eps << "; N=" << iterationsMax << std::endl;

    int i = 0;
    for (i = 0; i < iterationsMax; i++) {
        if (info) {
            std::cout << "Solver_BiCGSTAB: " << i << ": res = " << res << " tol = " << res / initres << std::endl;
        }
        if (res < eps) break;
        if (res > initres/eps) return -1;

        if (i == 0) Rhoi_1 = initres*initres;
        else Rhoi_1 = dotProduct(rr2, rr);
        if (fabs(Rhoi_1) < RhoMin) return -1;

        if (i == 0) pp = rr;
        else {
            betai_1 = (Rhoi_1*alphai_1) / (Rhoi_2*wi_1);
            linearCombination(pp, rr, betai_1, 1.0);
            // p = r + betai_1*(p - wi*v);
            linearCombination(pp, vv, 1.0, -wi_1*betai_1);
        }

        dd->matrixVectorProduct(pp, pp2);
        matrixA->matrixVectorProduct(pp2, vv);

        alphai = dotProduct(rr2, vv);
        if (fabs(alphai) < RhoMin) return -3;
        alphai = Rhoi_1 / alphai;

        ss = rr;
        // s = r - alphai*v;
        linearCombination(ss, vv, 1.0, -alphai);

        dd->matrixVectorProduct(ss, ss2);

        matrixA->matrixVectorProduct(ss2, tt);

        wi = dotProduct(tt, tt);
        if (fabs(wi) < RhoMin) return -4;
        wi = dotProduct(tt, ss) / wi;
        if (fabs(wi) < RhoMin) return -5;

        // x = x + alphai*p2 + wi*s2;
        linearCombination(solutionVector, pp2, 1.0, alphai);
        linearCombination(solutionVector, ss2, 1.0, wi);

        rr = ss;
        // r = s - wi*t;

        linearCombination(rr, tt, 1.0, -wi);

        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;

        res = sqrt(dotProduct(rr, rr));
    }
    if (info) {
        std::cout << "Solver_BiCGSTAB: outres: " << res << std::endl;
        std::cout << "Solver finished in " << i << " iterations, res = " << res << " tol=" << criterionTol << std::endl;
        return i;
    }
}

int solverBiCGSTAB_par(std::vector<double> & solutionVector, int sizeN, CompressedSparseRowMatrix * matrixA, std::vector<double> & vectorBB, double criterionTol, int iterationsMax) {
    DiagonalGenerator generatorDD;
    CompressedSparseRowMatrix * dd = generatorDD.generateDiagonalMatrix(matrixA);
    
    std::vector<double> pp(sizeN, 0), pp2(sizeN, 0), rr = vectorBB, rr2 = vectorBB, tt(sizeN, 0), vv(sizeN, 0), ss(sizeN, 0), ss2(sizeN, 0);

    double initres = sqrt(dotProduct_par(vectorBB, vectorBB)), res = initres;
    double mineps = 1E-15;
    double eps = std::max(criterionTol*initres, mineps);

    double Rhoi_1 = 1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, Rhoi_2 = 1.0, alphai_1 = 1.0, wi_1 = 1.0;
    double RhoMin = 1E-60;

    int info = 1;

    if (info) {
        std::cout << "Solver_BiCGSTAB: " << "initres: " << initres << "; eps: " << eps << "; N=" << iterationsMax << std::endl;
    }

    int i = 1;
    for (i = 0; i < iterationsMax; i++) {
        if (info) {
            std::cout << "Solver_BiCGSTAB: " << i << ": res = " << res << " tol = " << res / initres << std::endl;
        }
        if (res < eps) break;
        if (res > initres/eps) return -1;

        if (i == 0) Rhoi_1 = initres*initres;
        else Rhoi_1 = dotProduct_par(rr2, rr);
        if (fabs(Rhoi_1) < RhoMin) return -1;

        if (i == 0) pp = rr;
        else {
            betai_1 = (Rhoi_1*alphai_1) / (Rhoi_2*wi_1);
            linearCombination_par(pp, rr, betai_1, 1.0);
            // p = r + betai_1*(p - wi*v);
            linearCombination_par(pp, vv, 1.0, -wi_1*betai_1);
        }

        dd->matrixVectorProduct_par(pp, pp2);
        matrixA->matrixVectorProduct_par(pp2, vv);

        alphai = dotProduct_par(rr2, vv);
        if (fabs(alphai) < RhoMin) return -3;
        alphai = Rhoi_1 / alphai;

        ss = rr;
        // s = r - alphai*v;
        linearCombination_par(ss, vv, 1.0, -alphai);

        dd->matrixVectorProduct_par(ss, ss2);

        matrixA->matrixVectorProduct_par(ss2, tt);

        wi = dotProduct_par(tt, tt);
        if (fabs(wi) < RhoMin) return -4;
        wi = dotProduct_par(tt, ss) / wi;
        if (fabs(wi) < RhoMin) return -5;

        // x = x + alphai*p2 + wi*s2;
        linearCombination_par(solutionVector, pp2, 1.0, alphai);
        linearCombination_par(solutionVector, ss2, 1.0, wi);

        rr = ss;
        // r = s - wi*t;

        linearCombination_par(rr, tt, 1.0, -wi);

        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;

        res = sqrt(dotProduct_par(rr, rr));
    }
    if (info) {
        std::cout << "Solver_BiCGSTAB: outres: " << res << std::endl;
        std::cout << "Solver finished in " << i << " iterations, res = " << res << " tol=" << criterionTol << std::endl;
        return i;
    }
}

int main(int argc, char **argv)
{
    int nx = 1, ny = 1, nz = 1; // grid topology size for matrix generation
    double tol = 1e-06; // residual relative to the vector b norm
    int maxit = 100; // maximum number of iterations
    int nt = 1; // number of threads
    int qt = 0; // basic operation test flag

    parseArgs(argc, argv, nx, ny, nz, tol, maxit, nt, qt);
    //std::cout << nx << " " << ny << " " << nz << " " << tol << " " << maxit << " " << nt << " " << qt << std::endl;

    std::cout << "Testing BiCGSTAB solver for a 3D grid domain" << std::endl;
    std::cout << "N = " << nx*ny*nz << " (Nx=" << nx << ", Ny=" << ny << ", Nz=" << nz << ")" << std::endl;
    std::cout << "Aij = sin((double)i+j+1), i!=j" << std::endl;
    std::cout << "Aii = 1.1*sum(fabs(Aij))" << std::endl;
    std::cout << "Bi = sin((double)i)" << std::endl;
    std::cout << "tol = " << tol << std::endl << std::endl;

    if (qt) {
        std::cout << "This option switches on testing the speedup of basic operations." << std::endl;
        std::cout << "Run ./test to test the correctness of basic operations and its parallel versions." << std::endl << std::endl;

        testBasicOperationsSpeedup(nx, ny, nz);
    }

    RegularGridGenerator generatorMatrixA(nx, ny, nz);
    CompressedSparseRowMatrix * matrixA = generatorMatrixA.generateDoubleCSRMatrix();

    std::vector<double> vectorBB;
    for (int i = 0; i < matrixA->size(); i++) {
        vectorBB.push_back(sin(i));
    }    

    std::vector<double> solutionVector(matrixA->size(), 0);

    double timeSeqSolver = omp_get_wtime();
    std::cout << "Testing sequential solver:" << std::endl;
    solverBiCGSTAB(solutionVector, matrixA->size(), matrixA, vectorBB, tol, maxit);
    std::cout << "Sequential solver took: " << omp_get_wtime()-timeSeqSolver << std::endl;
    std::cout << std::endl;
    
    omp_set_dynamic(0);
    omp_set_num_threads(nt);
    
    double timeParSolver = omp_get_wtime();
    std::cout << "Testing parallel solver:" << std::endl;
    solverBiCGSTAB_par(solutionVector, matrixA->size(), matrixA, vectorBB, tol, maxit);
    std::cout << "Parallel solver took: " << omp_get_wtime()-timeParSolver << std::endl;
    
    // prints solution vector
    /*for (int i = 0; i < solutionVector.size(); i++) {
        std::cout.precision(4);
        std::cout << solutionVector[i] << " ";
    }
    std::cout << std::endl;
    */

    return 0;
}
