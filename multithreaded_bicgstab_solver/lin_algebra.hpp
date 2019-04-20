#include <vector>
#include <omp.h>

double dotProduct(std::vector<double> & x, std::vector<double> & y) {
    double dot = 0;
    for (int i = 0; i < x.size(); i++) {
        dot += x[i] * y[i];
    }
    return dot;
}

double dotProduct_par(std::vector<double> & x, std::vector<double> & y) {
    double dot = 0;
    #pragma omp parallel for reduction(+:dot)
    for (int i = 0; i < x.size(); i++) {
        dot += x[i] * y[i];
    }
    return dot;
}

std::vector<double> & linearCombination(std::vector<double> & x, std::vector<double> & y, double a, double b) {
    for (int i = 0; i < x.size(); i++) {
        x[i] = a*x[i] + b*y[i];
    }
    return x;
}

std::vector<double> & linearCombination_par(std::vector<double> & x, std::vector<double> & y, double a, double b) {
    #pragma omp parallel for
    for (int i = 0; i < x.size(); i++) {
        x[i] = a*x[i] + b*y[i];
    }
    return x;
}