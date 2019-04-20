#include <vector>
// #include <omp.h>
#include <mpi.h>

double MPI_dotProduct(std::vector<double> & x, std::vector<double> & y, int sizeRows, int sizeHalo) {
    double dotPart = 0, dot = 0;

    // #pragma omp parallel for reduction(+:dot)
    for (int i = 0; i < sizeRows; i++) {
        dotPart += x[i] * y[i];
    }

    MPI_Allreduce(&dotPart, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return dot;
}

std::vector<double> & MPI_linearCombination(std::vector<double> & x, std::vector<double> & y, double a, double b,
                                            int sizeRows, int sizeHalo) {
    // #pragma omp parallel for
    for (int i = 0; i < sizeRows; i++) {
        x[i] = a*x[i] + b*y[i];
    }
    return x;
}