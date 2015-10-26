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
 * @brief TwoSampleTTest
 * @details 
 *         See PermTesting.h for information
 */
arma::mat TwoSampleTTest(arma::mat group1, arma::mat group2)
{
    int nGroup1 = group1.n_rows;
    int nGroup2 = group2.n_rows;

    arma::mat group1_mean = arma::mean(group1); // 1 x Voxels vector
    arma::mat group2_mean = arma::mean(group2); // 1 x Voxels vector
    arma::mat group1_var = arma::var(group1); // 1 x Voxels vector
    arma::mat group2_var = arma::var(group2); // 1 x Voxels vector
    arma::mat mean_difference = group1_mean - group2_mean;
    arma::mat denominator = sqrt( (group1_var / (nGroup1-1)) + (group2_var / (nGroup2-1)) );

    arma::mat tStat = mean_difference / denominator;

    return tStat;
}

/**
 * @brief PermTestingIterative
 * @details 
 *         See PermTesting.h for information
 */
arma::mat PermTestingIterative(arma::mat data, 
                               int nPermutations, 
                               int nGroup1)
{
    std::string path = "/Users/felipegb94/PermTest/data/raw_adrc/indexMatrix.arma";
    int N; // Total number of subjects
    int nGroup2; // Total number of subjects in group 2
    int V; // Total number of statistics/voxels to be tested
    arma::mat indexMatrix;
    arma::mat maxT; // Maximum null distribution
    arma::mat group1; // Temporary vector 
    arma::mat group2; // Temporary vector 
    arma::mat tStat; // Temporary vector 

    /* Get dimensions */  
    N = data.n_rows; 
    V = data.n_cols;
    nGroup2 = N - nGroup1;
    indexMatrix = GetIndexMatrix(nPermutations,N);
    maxT = arma::zeros(nPermutations, 1);
    group1 = arma::zeros(nGroup1, V);
    group2 = arma::zeros(nGroup2, V);
    tStat = arma::zeros(1, V);

    arma::urowvec label_j(1,N); // Temporary vector 

    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): " << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;

    /* Permutation loop */
    #if OPENMP_ENABLED
        #pragma omp parallel for
    #endif
    for(int i = 0;i < nPermutations ;i++ )
    {
        std::cout << "Permutation " << i << std::endl;
        for(int j = 0; j < N; j++)
        {
            label_j(j) = indexMatrix(i, j);
        }
        group1 = data.rows(label_j(arma::span(0, nGroup1-1)));
        group2 = data.rows(label_j(arma::span(nGroup1, N-1)));   
        tStat = TwoSampleTTest(group1, group2);
        maxT(i, arma::span::all) = arma::max(tStat,1);
    }
    std::string prefix = "MaxT_IterativeCPU";
    SaveMaxT(maxT, nPermutations, prefix);


    return maxT;
}

/**
 * @brief PermTestingMatrix
 * @details 
 *         See PermTesting.h for information
 */
arma::mat PermTestingMatrix(arma::mat data, 
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
    int intervalSize; // Number of permutations done at once
    int numPasses; // nPermutations/intervalSize
    int start, end; // Start and end of current interval
    arma::cube permutationMatrices(nPermutations, N, 2);
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
    intervalSize = GetIntervalDimension(V, maxMemory);
    numPasses = nPermutations / intervalSize;
    g1Mean = arma::zeros(intervalSize, N);
    g2Mean = arma::zeros(intervalSize, N);
    g1Var = arma::zeros(intervalSize, N);
    g2Var = arma::zeros(intervalSize, N);
    tStatMatrix = arma::zeros(intervalSize, N);
    maxT = arma::zeros(nPermutations,1);


    std::cout << "Number of subjects (rows in data matrix): " << N << std::endl;
    std::cout << "Number of voxels per subject (cols in data matrix and cols in indexMatrix): ";
    std::cout << V << std::endl;
    std::cout << "Number of Permutations (rows in permutations matrix):" << nPermutations << std::endl;
    std::cout << "Size of group1 = " << nGroup1 << std::endl;
    std::cout << "Size of group2 = " << nGroup2 << std::endl;
    std::cout << "Rows in PermutationMatrices = " << permutationMatrices.n_rows << std::endl;
    std::cout << "Cols in PermutationMatrices = " << permutationMatrices.n_cols << std::endl;
    std::cout << "Interval Size = " << intervalSize << std::endl;
    std::cout << "Number of Passes = " << numPasses << std::endl;

    if((nPermutations % intervalSize) != 0)
    {
        fprintf(stderr, "Error: Wrong intervalSize \n");
        return maxT;
    }

    /* Permutation loop */
    #if OPENMP_ENABLED
        #pragma omp parallel for
    #endif
    for(int i = 0;i < numPasses;i++)
    {
        std::cout << "Curr Pass = " << i << " out of " << numPasses << std::endl;
        start = intervalSize * i;
        end = (intervalSize * i) + intervalSize - 1;

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




