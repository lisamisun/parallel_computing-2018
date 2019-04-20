#include <stdlib.h>
#include <stdio.h>

#include "parse_args.hpp"
#include "lin_algebra.hpp"
#include "generator.hpp"
#include "test_basic_operations.hpp"

int MPI_solverBiCGSTAB(std::vector<double> & solutionVector, int sizeRows, int sizeHalo, 
                       MPI_CompressedSparseRowMatrix * matrixA, std::vector<double> & vectorBB, 
                       double criterionTol, int iterationsMax, int procN) {
    MPI_DiagonalGenerator generatorDD;
    MPI_CompressedSparseRowMatrix * dd = generatorDD.MPI_generateDiagonalMatrix(matrixA);
    
    std::vector<double> pp(sizeRows+sizeHalo, 0), 
                        pp2(sizeRows+sizeHalo, 0), 
                        rr = vectorBB, 
                        rr2 = vectorBB, 
                        tt(sizeRows+sizeHalo, 0), 
                        vv(sizeRows+sizeHalo, 0), 
                        ss(sizeRows+sizeHalo, 0), 
                        ss2(sizeRows+sizeHalo, 0);

    double initres = sqrt(MPI_dotProduct(vectorBB, vectorBB, sizeRows, sizeHalo)), res = initres;
    double mineps = 1E-15;
    double eps = std::max(criterionTol*initres, mineps);

    double Rhoi_1 = 1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, Rhoi_2 = 1.0, alphai_1 = 1.0, wi_1 = 1.0;
    double RhoMin = 1E-60;

    int info = 1;

    if (procN == 0) {
        std::cout << "Solver_BiCGSTAB: " << ": initres: " << initres << "; eps: " << eps << "; N=" << iterationsMax << std::endl;
    }

    int i = 0;
    for (i = 0; i < iterationsMax; i++) {
        if ((info) && (procN == 0)) {
            std::cout << "Solver_BiCGSTAB: " << i << ": res = " << res << " tol = " << res / initres << std::endl;
        }
        if (res < eps) break;
        if (res > initres/eps) return -1;

        if (i == 0) Rhoi_1 = initres*initres;
        else Rhoi_1 = MPI_dotProduct(rr2, rr, sizeRows, sizeHalo);
        if (fabs(Rhoi_1) < RhoMin) return -1;

        if (i == 0) pp = rr;
        else {
            betai_1 = (Rhoi_1*alphai_1) / (Rhoi_2*wi_1);
            MPI_linearCombination(pp, rr, betai_1, 1.0, sizeRows, sizeHalo);
            // p = r + betai_1*(p - wi*v);
            MPI_linearCombination(pp, vv, 1.0, -wi_1*betai_1, sizeRows, sizeHalo);
        }

        dd->MPI_matrixVectorProduct(pp, pp2);

        MPI_update(pp2, matrixA->getRecieveRows(), matrixA->getSendRows(), matrixA->getGlobToLoc());
        matrixA->MPI_matrixVectorProduct(pp2, vv);

        alphai = MPI_dotProduct(rr2, vv, sizeRows, sizeHalo);
        if (fabs(alphai) < RhoMin) return -3;
        alphai = Rhoi_1 / alphai;

        ss = rr;
        // s = r - alphai*v;
        MPI_linearCombination(ss, vv, 1.0, -alphai, sizeRows, sizeHalo);

        dd->MPI_matrixVectorProduct(ss, ss2);

        MPI_update(ss2, matrixA->getRecieveRows(), matrixA->getSendRows(), matrixA->getGlobToLoc());
        matrixA->MPI_matrixVectorProduct(ss2, tt);

        wi = MPI_dotProduct(tt, tt, sizeRows, sizeHalo);
        if (fabs(wi) < RhoMin) return -4;
        wi = MPI_dotProduct(tt, ss, sizeRows, sizeHalo) / wi;
        if (fabs(wi) < RhoMin) return -5;

        // x = x + alphai*p2 + wi*s2;
        MPI_linearCombination(solutionVector, pp2, 1.0, alphai, sizeRows, sizeHalo);
        MPI_linearCombination(solutionVector, ss2, 1.0, wi, sizeRows, sizeHalo);

        rr = ss;
        // r = s - wi*t;

        MPI_linearCombination(rr, tt, 1.0, -wi, sizeRows, sizeHalo);

        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;

        res = sqrt(MPI_dotProduct(rr, rr, sizeRows, sizeHalo));
    }
    if ((info) && (procN == 0)) {
        std::cout << "Solver_BiCGSTAB: outres: " << res << std::endl;
        std::cout << "Solver finished in " << i << " iterations, res = " << res << " tol=" << criterionTol << std::endl;
    }
    return i;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int nProcs, procN;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // total number of processes in one switch
    MPI_Comm_rank(MPI_COMM_WORLD, &procN); // process number

    // parses arguments
    int nx = 1, ny = 1, nz = 1; // grid topology size for matrix generation
    double tol = 1e-06; // residual relative to the vector b norm
    int maxit = 100; // maximum number of iterations
    int nt = 1; // number of threads
    int qt = 0; // basic operation test flag
    int px = 1, py = 1, pz = 1; // number of parts for dividing the grid for each coordinate

    parseArgs(argc, argv, nx, ny, nz, tol, maxit, nt, qt, px, py, pz);
    //std::cout << nx << " " << ny << " " << nz << " " << tol << " " << maxit << " " << nt << " " << qt << std::endl;

    if (procN == 0) {
        std::cout << "Testing BiCGSTAB solver for a 3D grid domain" << std::endl;
        std::cout << "N = " << nx*ny*nz << " (Nx=" << nx << ", Ny=" << ny << ", Nz=" << nz << ")" << std::endl;
        std::cout << "P = " << px*py*pz << " (Px=" << px << ", Py=" << py << ", Pz=" << pz << ")" << std::endl;
        std::cout << "Aij = sin((double)i+j+1), i!=j" << std::endl;
        std::cout << "Aii = 1.1*sum(fabs(Aij))" << std::endl;
        std::cout << "Bi = sin((double)i)" << std::endl;
        std::cout << "tol = " << tol << std::endl << std::endl;
    }

    if (qt) {
        if (procN == 0) {
            std::cout << "Option 'qt' switches on testing the speedup of basic operations." << std::endl << std::endl;
        }

        MPI_testBasicOperationsTime(nx, ny, nz, px, py, pz, procN);
    }

    if (px*py*pz != nProcs) {
        if (procN == 0) {
            std::cout << "Number of parts for dividing the grid for each coordinate and number of MPI processe don't match!" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    /* omp_set_dynamic(0);
    omp_set_num_threads(nt);
    */

    MPI_RegularGridGenerator MPI_generatorMatrixA(nx, ny, nz, px, py, pz, procN);
    MPI_CompressedSparseRowMatrix * MPI_matrixA = MPI_generatorMatrixA.MPI_generateDoubleCSRMatrix();

    std::vector<double> MPI_vectorBB(MPI_matrixA->rowsHaloSize(), 0);
    for (int i = 0; i < MPI_matrixA->rowsSize(); i++) {
        MPI_vectorBB[i] = sin(MPI_matrixA->getGlobRowNumber(i));
    }    

    std::vector<double> MPI_solutionVector(MPI_matrixA->rowsHaloSize(), 0);

    MPI_Barrier(MPI_COMM_WORLD); // waiting for the initialization of all processes

    if (procN == 0) {
        std::cout << "Testing MPI parallel solver (without OpenMP):" << std::endl;
    }
    double procTimeMPISolver = MPI_Wtime();
    MPI_solverBiCGSTAB(MPI_solutionVector, MPI_matrixA->rowsSize(), MPI_matrixA->haloSize(), MPI_matrixA, MPI_vectorBB, tol, maxit, procN);
    procTimeMPISolver = MPI_Wtime() - procTimeMPISolver;
    double timeMPISolver = 0;
    MPI_Reduce(&procTimeMPISolver, &timeMPISolver, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (procN == 0) {
        std::cout << "MPI parallel solver took: " << timeMPISolver << "s" << std::endl;
    }

    // prints solution vector
    /*for (int i = 0; i < solutionVector.size(); i++) {
        std::cout.precision(4);
        std::cout << solutionVector[i] << " ";
    }
    std::cout << std::endl;
    */

    MPI_Finalize();
    return 0;
}
