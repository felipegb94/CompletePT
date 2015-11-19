#define ARMA_64BIT_WORD
#include "armadillo"
#include "PermTestingCPU.h"
#include "Utils/TimerCPU.h"



int main()
{
    std::string dataArmaPath = "/Users/felipegb94/PermTest/data/face/Data_face.arma";

    int N;
    int V;
    int nPermutations = 4096;
    int i;
    double maxMemory = 5e8; // Max number of bytes to use

    arma::mat data;
    data.load(dataArmaPath);

    arma::mat maxT = OneSamplePermTestingCPU(data, nPermutations, maxMemory);



    return 0;
}