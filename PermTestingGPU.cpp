#include <iostream>
#include <vector>
#include <stdlib.h>

#include "arrayfire.h"
#include "armadillo"
#include "Utils/PermTestingUtils.h"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"

float* PermTestingGPU(float* dataHost, float* permutationMatrix1Host, float* permutationMatrix2Host,
                      int N, int V, int nPermutations, int nGroup1, double maxMemory)
{
    std::string path1 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix1.arma";
    std::string path2 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix2.arma";
    std::cout << "Starting Permutation Testing with ArrayFire!" << std::endl;
    std::cout << "Environment Information: " << std::endl;
    af::info();

    arma::mat permutationMatrices = GetPermutationMatrices(nPermutations, N, nGroup1);
    arma::mat permutationMatrix1Host = permutationMatrices.slice(0);
    arma::mat permutationMatrix2Host = permutationMatrices.slice(1);

    af::array dataDevice(V,N,dataHost);
    af::array dataSquaredDevice = dataDevice * dataDevice;
    af::array permutationMatrix1Device(N,nPermutations,permutationMatrix1Host);
    af::array permutationMatrix2Device(N,nPermutations,permutationMatrix2Host);

    int intervalSize = GetIntervalDimension(V, maxMemory);
    int numPasses = nPermutations/intervalSize;

    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Size of group1 = " << nGroup1 << std::endl;
    std::cout << "Size of group2 = " << nGroup2 << std::endl;
    std::cout << "Rows in PermutationMatrices = " << N << std::endl;
    std::cout << "Cols in PermutationMatrices = " << nPermutations << std::endl;
    std::cout << "Interval Size = " << intervalSize << std::endl;
    std::cout << "Number of Passes = " << numPasses << std::endl;


    af::array g1Mean(N, intervalSize);
    af::array g2Mean(N, intervalSize);
    af::array g1Var(N, intervalSize);
    af::array g2Var(N, intervalSize);
    af::array tStatMatrix(N, intervalSize);
    af::array MaxTDevice(1,nPermutations);

    int start, end;
    for(int i = 0;i < numPasses;i++)
    {
        start = intervalSize * i;
        end = (intervalSize * i) + intervalSize - 1;
        std::cout << "Curr Pass = " << i << std::endl;
        g1Mean = af::matmul(dataDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / N_g1 ;
        g2Mean = af::matmul(dataDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / N_g2;
        g1Var = (af::matmul(dataSquaredDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / (N_g1)) - (g1Mean*g1Mean);

        g2Var = (af::matmul(dataSquaredDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / (N_g2)) - (g2Mean*g2Mean); 

        tStatMatrix = (g1Mean - g2Mean) / (af::sqrt((g1Var/(N_g1-1)) + (g2Var/(N_g2-1))));

        MaxTDevice(af::seq(start,end)) = af::max(tStatMatrix,0);

    }
    
    float * MaxTHost = MaxTDevice.host<float>();
    arma::mat MaxT(1,nPermutations);

    for(int i = 0;i < nPermutations;i++)
    {
        MaxT(i) = MaxTHost[i];
    }
    MaxT.save("MaxT_ArrayFire_10000.arma");
    MaxT.save("MaxT_ArrayFire_10000.ascii", arma::raw_ascii);




    return MaxTHost;

}