#ifndef PERMTESTINGSHARED_H
#define PERMTESTINGSHARED_H
#include "armadillo"

/**
 * @brief TwoSampleGetPermutationMatrices
 * @details 
 *         Generate a matrix of 1's and 0's corresponding to group1. In each row, if an element is 
 *         a 1 it means that the row in the data matrix at that same index makes part of group1. 
 *         If it is 0 it means it belongs to group2. So each row has a different combination of 
 *         1's and 0's.
 * @param nPermutations Number of permutations
 * @param N Total number of indeces to permute
 * 
 * @return indexMatrix - nPermutationsxN matrix
 */
arma::cube TwoSampleGetPermutationMatrices(int nPermutations, int N, int nGroup1);

     
/**
 * @brief OneSampleGetPermutationMatrix
 * @details 
 */
arma::mat OneSampleGetPermutationMatrix(int nPermutations, int N);

#endif