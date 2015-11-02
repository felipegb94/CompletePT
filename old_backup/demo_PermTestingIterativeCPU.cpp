#include <iostream>
#define ARMA_64BIT_WORD
#include "armadillo"

#include "PermTestingCPU.h"
#include "Utils/TimerCPU.h"


int main()
{
    std::string dataArmaPath = "/Users/felipegb94/PermTest/data/raw_adrc/adrc_raw.arma";
    int nGroup1 = 25;
    int nPermutations = 100;

    int N;
    int V;
    arma::mat data;
    arma::mat indexMatrix;
    arma::mat T;
    arma::mat MaxT;

    data.load(dataArmaPath);
    N = data.n_rows; 
    V = data.n_cols;

    std::cout << "Number of subjects (rows in data matrix): ";
    std::cout << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix): ";
    std::cout << nPermutations << std::endl;
    std::cout << "Rows in IndexMatrix = " << indexMatrix.n_rows << std::endl;
    std::cout << "Cols in IndexMatrix = " << indexMatrix.n_cols << std::endl;

    CpuTimer t;
    t.Start();
    t.Start_cputimer();
    arma::mat maxT = PermTestingIterative(data, nPermutations, nGroup1);
    t.Stop();
    t.Stop_cputimer();

    t.print();
    return 0;
}



