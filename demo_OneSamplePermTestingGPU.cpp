#include <string>
#include "armadillo"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"
#include "PermTestingGPU.h"


int main()
{
    std::string dataArmaPath = "/Users/felipegb94/PermTest/data/face/Data_face.arma";
    arma::mat data;
    data.load(dataArmaPath);
    float * dataHost = ArmaToArray(data);
    int N = data.n_rows;
    int V = data.n_cols;
    int nPermutations = 4096;
    int i;
    double maxMemory = 5e8; // Max number of bytes to use


    float* maxT = OneSamplePermTestingGPU(dataHost, nPermutations, N, V, maxMemory);



    return 0;
}