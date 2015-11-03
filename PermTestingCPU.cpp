#include <math.h>
#include "sstream"
#include "vector"

#include "armadillo"
#include "Utils/PermTestingUtils.h"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"

#if OPENMP_ENABLED
    #include <omp.h>
#endif


/**
 * @brief PermTestingCPU
 * @details 
 *         See PermTesting.h for information
 */
arma::mat PermTestingCPU(arma::mat data, 
                            int nPermutations,
                            int nGroup1, 
                            double maxMemory)
{
    std::string path1 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix1.arma";
    std::string path2 = "/Users/felipegb94/PermTest/data/raw_adrc/permutationMatrix2.arma";
    std::cout << "Starting PermTestingMatrix..." << std::endl;
    int N; // Total number of subjects
    int nGroup2; // Total number of subjects in group 2
    int V; // Total number of statistics/voxels to be tested
    int nPermutationsPerIteration; // Number of permutations done at once
    int numIterations; // nPermutations/nPermutationsPerIteration
    int lastIteration;
    int start, end; // Start and end of current interval

    arma::cube permutationMatrices;
    arma::mat permutationMatrix1;
    arma::mat permutationMatrix2; 
    arma::mat maxT; // Maximum null distribution
    arma::mat dataSquared;
    arma::mat g1Mean;
    arma::mat g2Mean;
    arma::mat g1Var;
    arma::mat g2Var;
    arma::mat tStatMatrix;

    /* Set constant values and allocate memory */   
    N = data.n_rows; 
    V = data.n_cols;
    nGroup2 = N - nGroup1;
    permutationMatrices = GetPermutationMatrices(nPermutations, N, nGroup1);
    permutationMatrix1 = permutationMatrices.slice(0);
    permutationMatrix2 = permutationMatrices.slice(1);
    dataSquared = data % data;
    nPermutationsPerIteration = GetIntervalDimension(V, maxMemory);
    lastIteration = nPermutations % nPermutationsPerIteration;
    numIterations = floor(nPermutations/nPermutationsPerIteration);
    g1Mean = arma::zeros(nPermutationsPerIteration, N);
    g2Mean = arma::zeros(nPermutationsPerIteration, N);
    g1Var = arma::zeros(nPermutationsPerIteration, N);
    g2Var = arma::zeros(nPermutationsPerIteration, N);
    tStatMatrix = arma::zeros(nPermutationsPerIteration, N);
    maxT = arma::zeros(nPermutations,1);


    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Size of group1 = " << nGroup1 << std::endl;
    std::cout << "Size of group2 = " << nGroup2 << std::endl;
    std::cout << "Rows in PermutationMatrices = " << permutationMatrices.n_rows << std::endl;
    std::cout << "Cols in PermutationMatrices = " << permutationMatrices.n_cols << std::endl;
    std::cout << "Interval Size = " << nPermutationsPerIteration << std::endl;
    std::cout << "Number of Passes = " << numIterations << std::endl;


    int i = 0;
    /* Permutation loop */
    #if OPENMP_ENABLED
        #pragma omp parallel for
    #endif
    for(i = 0;i < numIterations;i++)
    {
        start = nPermutationsPerIteration * i;
        end = (nPermutationsPerIteration * i) + nPermutationsPerIteration - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);

        g1Mean = (permutationMatrix1(arma::span(start,end), arma::span::all) * data) / nGroup1;
        g2Mean = (permutationMatrix2(arma::span(start,end), arma::span::all) * data) / nGroup2;
        g1Var = ((permutationMatrix1(arma::span(start,end), arma::span::all) * dataSquared) / (nGroup1)) - (g1Mean % g1Mean); 
        g2Var = ((permutationMatrix2(arma::span(start,end), arma::span::all) * dataSquared) / (nGroup2)) - (g2Mean % g2Mean); 

        tStatMatrix = (g1Mean - g2Mean) / (sqrt((g1Var/(nGroup1-1)) + (g2Var/(nGroup2-1))));
        maxT(arma::span(start,end),arma::span::all) = arma::max(tStatMatrix,1);
    }
    if(lastIteration != 0)
    {
        start = nPermutationsPerIteration * i;
        end = nPermutations - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);

        g1Mean = (permutationMatrix1(arma::span(start,end), arma::span::all) * data) / nGroup1;
        g2Mean = (permutationMatrix2(arma::span(start,end), arma::span::all) * data) / nGroup2;
        g1Var = ((permutationMatrix1(arma::span(start,end), arma::span::all) * dataSquared) / (nGroup1)) - (g1Mean % g1Mean); 
        g2Var = ((permutationMatrix2(arma::span(start,end), arma::span::all) * dataSquared) / (nGroup2)) - (g2Mean % g2Mean); 

        tStatMatrix = (g1Mean - g2Mean) / (sqrt((g1Var/(nGroup1-1)) + (g2Var/(nGroup2-1))));
        maxT(arma::span(start,end),arma::span::all) = arma::max(tStatMatrix,1);   
    }
    std::string prefix = "MaxT_CPU";
    SaveMaxT(maxT, nPermutations, prefix);

    return maxT;
}

arma::mat OneSampleGetPermutationMatrix(int nPermutations, int N)
{
    std::cout << "Getting one sample permMatrix " << std::endl;
    arma::arma_rng::set_seed_random();  // set the seed to a random value

    int cutoff;
    bool cutoffOne = false;
    arma::mat indexList; 
    arma::mat permutationMatrix(nPermutations, N, arma::fill::ones);

    indexList = arma::linspace<arma::mat>(0, N-1, N);
    indexList = arma::shuffle(indexList);
    for(int i = 0;i < nPermutations;i++)
    {
        cutoff = indexList(0,0);
        indexList = arma::shuffle(indexList);
        std::cout << "cutoff = " << cutoff <<std::endl;
        indexList.print();
        for(int j = 0;j < cutoff;j++)
        {
            permutationMatrix(i, indexList(j,0)) = -1;
        }
        permutationMatrix(i,arma::span::all).print();
    }

    return permutationMatrix;
}

arma::mat OneSamplePermTestingCPU(arma::mat data, 
                                  int nPermutations,
                                  double maxMemory);
{
    int N; // Total number of subjects
    int V; // Total number of statistics/voxels to be tested
    int nPermutationsPerIteration; // Number of permutations done at once
    int numIterations; // nPermutations/nPermutationsPerIteration
    int lastIteration;
    int start, end; // Start and end of current interval

    arma::mat permutationMatrix;
    arma::mat maxT; // Maximum null distribution
    arma::mat dataSquared;
    arma::mat gMean;
    arma::mat gVar;
    arma::mat tStatMatrix;

    /* Set constant values and allocate memory */   
    N = data.n_rows; 
    V = data.n_cols;
    permutationMatrices = OneSampleGetPermutationMatrices(nPermutations, N);
    dataSquared = data % data;
    nPermutationsPerIteration = GetIntervalDimension(V, maxMemory);
    lastIteration = nPermutations % nPermutationsPerIteration;
    numIterations = floor(nPermutations/nPermutationsPerIteration);
    gMean = arma::zeros(nPermutationsPerIteration, N);
    gVar = arma::zeros(nPermutationsPerIteration, N);
    tStatMatrix = arma::zeros(nPermutationsPerIteration, N);
    maxT = arma::zeros(nPermutations,1);

    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Rows in PermutationMatrices = " << permutationMatrices.n_rows << std::endl;
    std::cout << "Cols in PermutationMatrices = " << permutationMatrices.n_cols << std::endl;
    std::cout << "Interval Size = " << nPermutationsPerIteration << std::endl;
    std::cout << "Number of Passes = " << numIterations << std::endl;

    return maxT;
}



