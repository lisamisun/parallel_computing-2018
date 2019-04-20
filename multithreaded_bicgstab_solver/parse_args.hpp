#include <stdlib.h>
#include <cstdlib>
#include <string>

void parseArgs(int argc, char **argv, int & nx, int & ny, int & nz, double & tol, int & maxit, int & nt, int & qt) {
    for (int i = 1; i < argc; i++) {
        std::string argument(argv[i]);
        if (argument.find("nx=") != std::string::npos) {
            nx = std::atoi(argument.substr(3).c_str());
        }
        if (argument.find("ny=") != std::string::npos) {
            ny = std::atoi(argument.substr(3).c_str());
        }
        if (argument.find("nz=") != std::string::npos) {
            nz = std::atoi(argument.substr(3).c_str());
        }
        if (argument.find("tol=") != std::string::npos) {
            tol = std::atof(argument.substr(4).c_str());
        }
        if (argument.find("maxit=") != std::string::npos) {
            maxit = std::atoi(argument.substr(6).c_str());
        }
        if (argument.find("nt=") != std::string::npos) {
            nt = std::atoi(argument.substr(3).c_str());
        }
        if (argument.find("qt") != std::string::npos) {
            qt = 1;
        }
    }
}
