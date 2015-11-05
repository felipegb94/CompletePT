#ifndef PERMTESTINGARMADILLOUTILS_H
#define PERMTESTINGARMADILLOUTILS_H

#include <string>
#include "armadillo"

/**
 * @brief SaveMaxT
 * @details 
 *         Saves maxT in arma and ascii format. If matrix == 1, it means it is saving a maxT that 
 *         was output from PermTestingMatrix. Else it was output from PermTestingIterative
 *         
 *         filename formats = "maxT_Matrix_numPermutations.arma"
 *                            "maxT_Matrix_numPermutations.ascii"
 *                            "maxT_Iterative_numPermutations.arma"
 *                            "maxT_Iterative_numPermutations.ascii"
 * 
 * @param maxT - maxT to save.
 * @param nPermutations - Used for the filename
 * @param matrix - Where did this maxT came from
 */
void SaveMaxT(arma::mat maxT, int nPermutations, std::string prefix);

/**
 * @brief ArmaToArray
 * @details 
 *         Convert an armadillo matrix into a float array.
 * return array
 */
float* ArmaToArray(arma::mat matrix);

/**
 * @brief GetIndexMatrix
 * @details 
 *         Generate a matrix where each row is a permutation of the list [0,1,2...N].
 *         This function is needed for PermTestingIterative.
 * @param nPermutations Number of permutations
 * @param N Total number of indeces to permute
 * 
 * @return indexMatrix - nPermutationsxN matrix
 */
arma::mat GetIndexMatrix(int nPermutations, int N);

/**
 * @brief LoadIndexMatrix
 * @details 
 *         LoadIndexMatrix from specified path
 * @param path - string containing the path to the matrix
 * 
 * @return indexMatrix - nPermutationsxN matrix
 */
arma::mat LoadIndexMatrix(std::string path, bool isArma);


/**
 * @brief LoadPermutationMatrices
 * @details 
 *         Load permutation matrix 1 and 2
 * @param path1 Path to permutationMatrix1
 * @param path2 Path to permutationMatrix2
 * 
 * @return indexMatrix - nPermutationsxN matrix
 */
arma::cube LoadPermutationMatrices(std::string path1, std::string path2, bool isArma);

#endif