#include <string>
#include "armadillo"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"
#include "PermTestingGPU.h"

int main()
{
    std::string dataArmaPath = "/home/felipe/PermTest/data/raw_adrc/adrc_raw.arma";
    arma::mat data;
    data.load(dataArmaPath);
    float * dataHost = ArmaToArray(data);

    int N = data.n_rows;
    int V = data.n_cols;
    int nGroup1 = 25;
    int nGroup2 = N - nGroup1;
    int nPermutations = 10000;
    int i;
    double maxMemory = 1e7; // Max number of bytes to use
    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Size of group1 = " << nGroup1 << std::endl;
    std::cout << "Size of group2 = " << nGroup2 << std::endl;

    TwoSamplePermTestingGPU(dataHost, nPermutations, N, V, nGroup1, maxMemory);

}