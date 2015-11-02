#include <iostream>
#include <vector>
#include <stdlib.h>

#include "arrayfire.h"
#include "armadillo"
#include "Utils/PermTestingUtils.h"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"

float* PermTestingGPU(float* dataHost, int nPermutations, int N, int V, int nGroup1, double maxMemory)
{
    std::string path1 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix1.arma";
    std::string path2 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix2.arma";
    std::cout << "Starting Permutation Testing with ArrayFire!" << std::endl;
    std::cout << "Environment Information: " << std::endl;
    af::info();

    int nGroup2 = N - nGroup1;
    arma::cube permutationMatrices = GetPermutationMatrices(nPermutations, N, nGroup1);
    arma::mat permutationMatrix1ArmaHost = permutationMatrices.slice(0);
    arma::mat permutationMatrix2ArmaHost = permutationMatrices.slice(1);
    float* permutationMatrix1Host = ArmaToArray(permutationMatrix1ArmaHost);
    float* permutationMatrix2Host = ArmaToArray(permutationMatrix2ArmaHost);

    int nPermutationsPerIteration = GetIntervalDimension(V, maxMemory);
    int lastIteration = nPermutations % nPermutationsPerIteration;
    int numIterations = floor(nPermutations/nPermutationsPerIteration);

    af::array dataDevice(V,N,dataHost);
    af::array dataSquaredDevice = dataDevice * dataDevice;
    af::array permutationMatrix1Device(N,nPermutations,permutationMatrix1Host);
    af::array permutationMatrix2Device(N,nPermutations,permutationMatrix2Host);


    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Size of group1 = " << nGroup1 << std::endl;
    std::cout << "Size of group2 = " << nGroup2 << std::endl;
    std::cout << "Rows in PermutationMatrices = " << N << std::endl;
    std::cout << "Cols in PermutationMatrices = " << nPermutations << std::endl;
    std::cout << "Interval Size = " << nPermutationsPerIteration << std::endl;
    std::cout << "Number of Passes = " << numIterations << std::endl;


    af::array g1Mean(N, nPermutationsPerIteration);
    af::array g2Mean(N, nPermutationsPerIteration);
    af::array g1Var(N, nPermutationsPerIteration);
    af::array g2Var(N, nPermutationsPerIteration);
    af::array tStatMatrix(N, nPermutationsPerIteration);
    af::array MaxTDevice(1,nPermutations);
    int i = 0;
    int start, end;
    for(i = 0;i < numIterations;i++)
    {
        start = nPermutationsPerIteration * i;
        end = (nPermutationsPerIteration * i) + nPermutationsPerIteration - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);

        g1Mean = af::matmul(dataDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / nGroup1 ;
        g2Mean = af::matmul(dataDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / nGroup2;
        g1Var = (af::matmul(dataSquaredDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / (nGroup1)) - (g1Mean*g1Mean);

        g2Var = (af::matmul(dataSquaredDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / (nGroup2)) - (g2Mean*g2Mean); 

        tStatMatrix = (g1Mean - g2Mean) / (af::sqrt((g1Var/(nGroup1-1)) + (g2Var/(nGroup2-1))));

        MaxTDevice(af::seq(start,end)) = af::max(tStatMatrix,0);

    }
    if(lastIteration != 0)
    {
        start = nPermutationsPerIteration * i;
        end = nPermutations - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);
        
        g1Mean = af::matmul(dataDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / nGroup1 ;
        g2Mean = af::matmul(dataDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / nGroup2;
        g1Var = (af::matmul(dataSquaredDevice, permutationMatrix1Device(af::span, af::seq(start, end))) / (nGroup1)) - (g1Mean*g1Mean);

        g2Var = (af::matmul(dataSquaredDevice, permutationMatrix2Device(af::span, af::seq(start, end))) / (nGroup2)) - (g2Mean*g2Mean); 

        tStatMatrix = (g1Mean - g2Mean) / (af::sqrt((g1Var/(nGroup1-1)) + (g2Var/(nGroup2-1))));

        MaxTDevice(af::seq(start,end)) = af::max(tStatMatrix,0);
    }

    float * MaxTHost = MaxTDevice.host<float>();
    arma::mat MaxT(1,nPermutations);
    for(int i = 0;i < nPermutations;i++)
    {
        MaxT(i) = MaxTHost[i];
    }
    std::string prefix = "MaxT_GPU";
    SaveMaxT(MaxT, nPermutations, prefix);

    return MaxTHost;
}