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