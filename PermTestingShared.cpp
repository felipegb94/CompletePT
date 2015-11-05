#include "armadillo"


/**
 * @brief GetPermutationMatrices
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
arma::cube TwoSampleGetPermutationMatrices(int nPermutations, int N, int nGroup1)
{
    arma::arma_rng::set_seed_random();  // set the seed to a random value

    arma::cube permutationMatrices(nPermutations, N, 2, arma::fill::zeros);

    arma::mat permutationMatrix1(nPermutations, N, arma::fill::zeros);

    arma::mat permutationMatrix2 = arma::ones(nPermutations, N);

    arma::mat indexList; 

    indexList = arma::linspace<arma::mat>(0, N-1, N);


    for(int i = 0;i < nPermutations;i++)
    {   
        indexList = arma::shuffle(indexList);
        for(int j = 0;j < nGroup1;j++)
        {
            permutationMatrix1(i,indexList(j)) = 1;
        }
    }
    permutationMatrix2 = permutationMatrix2 - permutationMatrix1;
    permutationMatrices.slice(0) = permutationMatrix1;
    permutationMatrices.slice(1) = permutationMatrix2;

    return permutationMatrices;
}

arma::mat OneSampleGetPermutationMatrix(int nPermutations, int N)
{
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
        for(int j = 0;j < cutoff;j++)
        {
            permutationMatrix(i, indexList(j,0)) = -1;
        }
    }

    return permutationMatrix;
}