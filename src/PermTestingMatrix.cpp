
#define ARMA_64BIT_WORD
#include "armadillo"

#include "PermTestingCPU.h"



int main()
{
    std::string dataArmaPath = "/Users/felipegb94/PermTest/data/raw_adrc/adrc_raw.arma";

    int nGroup1 = 25;
    int nGroup2;
    int N;
    int V;
    int nPermutations = 10000;
    int i;
    double maxMemory = 5e8; // Max number of bytes to use

    arma::mat data;
    data.load(dataArmaPath);

    PermTestingMatrix(data, nPermutations, nGroup1, maxMemory);

    return 0;
}