#include <math.h>
#include "sstream"
#include "vector"

#include "armadillo"
#include "Utils/PermTestingUtils.h"
#include "ArmadilloUtils/PermTestingArmadilloUtils.h"
#include "PermTestingShared.h"

#if OPENMP_ENABLED
    #include <omp.h>
#endif


/**
 * @brief PermTestingCPU
 * @details 
 *         See PermTesting.h for information
 */
arma::mat TwoSamplePermTestingCPU(arma::mat data, 
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
    std::cout << "testing GetPerm " <<std::endl;
    permutationMatrices = TwoSampleGetPermutationMatrices(nPermutations, N, nGroup1);
    std::cout << "testing GetPerm " <<std::endl;

    permutationMatrix1 = permutationMatrices.slice(0);
    permutationMatrix2 = permutationMatrices.slice(1);
    dataSquared = data % data;
    nPermutationsPerIteration = GetIntervalDimension(V, maxMemory);
    lastIteration = nPermutations % nPermutationsPerIteration;
    numIterations = floor(nPermutations/nPermutationsPerIteration);
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



arma::mat OneSamplePermTestingCPU(arma::mat data, 
                                  int nPermutations,
                                  double maxMemory)
{
    std::string permMatrixPath = "/Users/felipegb94/PermTest/data/face/PermutationMatrix.arma";
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
    //permutationMatrix.load(permMatrixPath);
    permutationMatrix = OneSampleGetPermutationMatrix(nPermutations, N);
    dataSquared = data % data;
    nPermutationsPerIteration = GetIntervalDimension(V, maxMemory);
    lastIteration = nPermutations % nPermutationsPerIteration;
    numIterations = floor(nPermutations/nPermutationsPerIteration);
    maxT = arma::zeros(nPermutations,1);

    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Rows in permutationMatrix = " << permutationMatrix.n_rows << std::endl;
    std::cout << "Cols in permutationMatrix = " << permutationMatrix.n_cols << std::endl;
    std::cout << "Interval Size = " << nPermutationsPerIteration << std::endl;
    std::cout << "Number of Passes = " << numIterations << std::endl;
    std::cout << "Last Iteration = " << lastIteration << std::endl;

    int i = 0;

    arma::mat varTerm1 = repmat(arma::sum(dataSquared,0), nPermutationsPerIteration,1)/N;
    arma::mat varTerm1LastItr = repmat(arma::sum(dataSquared,0), lastIteration, 1)/N;
    std::cout << "VarTerm rows " << varTerm1LastItr.n_rows << std::endl;
    std::cout << "VarTerm cols " << varTerm1LastItr.n_cols << std::endl;
    /* Permutation loop */
    #if OPENMP_ENABLED
        #pragma omp parallel for
    #endif
    for(i = 0;i < numIterations;i++)
    {
        start = nPermutationsPerIteration * i;
        end = (nPermutationsPerIteration * i) + nPermutationsPerIteration - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);

        gMean = (permutationMatrix(arma::span(start,end), arma::span::all) * data) / N;
        gVar = (varTerm1) - (gMean % gMean); 
        tStatMatrix = (gMean) / sqrt(gVar/(N-1));

        maxT(arma::span(start,end),arma::span::all) = arma::max(tStatMatrix,1);
    }
    if(lastIteration != 0)
    {
        start = nPermutationsPerIteration * i;
        end = nPermutations - 1;
        printf("Iteration %d , start %d, end %d of %d \n", i, start, end, nPermutations-1);

        gMean = (permutationMatrix(arma::span(start,end), arma::span::all) * data) / N;
        gVar = varTerm1LastItr - (gMean % gMean); 
        tStatMatrix = (gMean) / sqrt(gVar/(N-1));

        maxT(arma::span(start,end),arma::span::all) = arma::max(tStatMatrix,1);
    }

    std::string prefix = "OneSampleMaxT_CPU";
    SaveMaxT(maxT, nPermutations, prefix);

    return maxT;
}



