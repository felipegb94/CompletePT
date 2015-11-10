#include <vector>

#define ARMA_64BIT_WORD
#include <string>
#include "armadillo"
#include "PermTestingCPU.h"
#include "Utils/TimerCPU.h"



int main()
{
    std::string dataArmaPath = "/home/felipe/PermTest/data/ADRC/ADRC_100_50_50.arma";
    int numRuns = 8;
    int nPermutations[] = {2500, 5000, 10000, 20000, 40000, 80000, 160000, 320000};
    int nGroup1 = 50;
    int nGroup2;
    int N;
    int V;
    int i;
    double maxMemory = 5e8; // Max number of bytes to use

    arma::mat data;
    arma::mat timings = arma::zeros(2, numRuns);
    data.load(dataArmaPath);

    CpuTimer t;
    std::string prefix = "timing_CompletePT_ADRC";

    for(int i = 0;i < numRuns;i++)
    {
        t.Start();
        t.Start_cputimer();
        TwoSamplePermTestingCPU(data, nPermutations[i], nGroup1, maxMemory);
        t.Stop();
        t.Stop_cputimer();
        t.print();
        timings(0,i) = t.Elapsed();
        timings(1,i) = t.Elapsed_cputimer();

        std::stringstream armaFileName;
        std::stringstream asciiFileName;
        armaFileName << prefix << "_" << nPermutations[i] << ".arma";
        asciiFileName << prefix << "_" << nPermutations[i] << ".ascii";  
        timings.save(armaFileName.str());
        timings.save(asciiFileName.str(), arma::raw_ascii);
    }





    return 0;
}