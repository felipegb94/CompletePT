#include <iostream>
#include <math.h>
#include "sstream"
#include "vector"

#include "armadillo"


/**
 * @brief SaveMaxT
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
void SaveMaxT(arma::mat maxT, int nPermutations, std::string prefix)
{
    std::stringstream armaFileName;
    std::stringstream asciiFileName;

    armaFileName << prefix << "_" << nPermutations << ".arma";
    asciiFileName << prefix << "_" << nPermutations << ".ascii";
    
    maxT.save(armaFileName.str());
    maxT.save(asciiFileName.str(), arma::raw_ascii);
}

/**
 * @brief ArmaToArray
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
float* ArmaToArray(arma::mat matrix)
{
    int r = matrix.n_rows;
    int c = matrix.n_cols;
    float *arr = new float[r * c];

    for(int i = 0;i < r;i++)
    {
        for(int j = 0; j < c;j++)
        {
            arr[i*c + j] = matrix(i,j);
        }
    }     
    return arr;
}

/**
 * @brief GetIndexMatrix
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
arma::mat GetIndexMatrix(int nPermutations, int N)
{
    arma::arma_rng::set_seed_random();  // set the seed to a random value
    arma::mat indexMatrix;
    arma::mat indexList;

    indexList = arma::linspace<arma::mat>(0, N-1, N);
    indexMatrix = arma::zeros(N, nPermutations);

    for(int i = 0;i < nPermutations;i++)
    {
        indexList = arma::shuffle(indexList);
        indexMatrix(arma::span::all,i) = indexList;
    }

    return indexMatrix.t();
}

/**
 * @brief LoadIndexMatrix
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
arma::mat LoadIndexMatrix(std::string path, bool isArma)
{
    arma::mat indexMatrix;
    if(isArma)
    {
        indexMatrix.load(path);
    }
    else
    {
        indexMatrix.load(path, arma::raw_ascii);
    }
    return indexMatrix;
}

/**
 * @brief GetPermutationMatrices
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
arma::cube GetPermutationMatrices(int nPermutations, int N, int nGroup1)
{
    arma::arma_rng::set_seed_random();  // set the seed to a random value

    arma::cube permutationMatrices(nPermutations, N, 2, arma::fill::zeros);
    arma::mat permutationMatrix1(nPermutations, N, arma::fill::zeros);
    arma::mat permutationMatrix2(nPermutations, N, arma::fill::ones);
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

/**
 * @brief LoadPermutationMatrices
 * @details 
 *         See PermTestingArmadilloUtils.h for information
 */
arma::cube LoadPermutationMatrices(std::string path1, std::string path2, bool isArma)
{
    std::cout << "Loading PermutationMatrices..." << std::endl;
    arma::mat permutationMatrix1;
    arma::mat permutationMatrix2;

    if(isArma)
    {
        permutationMatrix1.load(path1);
        permutationMatrix2.load(path2);
    }
    else
    {
        permutationMatrix1.load(path1, arma::raw_ascii);
        permutationMatrix2.load(path2, arma::raw_ascii);
    }
    std::cout << "PM1 rows = " << permutationMatrix1.n_rows << std::endl;
    std::cout << "PM1 cols = " << permutationMatrix1.n_cols << std::endl;

    arma::cube permutationMatrices(permutationMatrix1.n_rows, permutationMatrix1.n_cols, 2);
    permutationMatrices.slice(0) = permutationMatrix1;
    permutationMatrices.slice(1) = permutationMatrix2;

    return permutationMatrices;   
}