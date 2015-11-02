#ifndef PERMTESTINGCPU_H
#define PERMTESTINGCPU_H
#include "armadillo"

/**
 * @brief TwoSampleTTest
 * @details 
 *         Performs two sample t-test between group 1 and group 2
 * 
 * @param group1 nGroup1 x V vector
 * @param group2 nGroup1 x V vector
 * 
 * @return tStat 1xV vector: Calculated t-stat for each voxel.
 */
arma::mat TwoSampleTTest(arma::mat group1, arma::mat group2);

/**
 * @brief PermTestingIterative
 * @details 
 *         The following method performs permutation testing to find group differences between 2 
 *         different groups using a two-sample t-test as the test statistic.
 * 
 * @param data NxV matrix: N = Total number of subjects. 
 *                         V = Total number of statistics/voxels per subject
 * @param indexMatrix nPermutationsxV Matrix: 
 *                         nPermutations = Total number of permutations to be done
 *                         V = Same as in data.
 * @param nGroup1 = Total number of subjects in group1. N-nGroup1 = nGroup2.
 * @return maxT numPermutationsx1 vector: Maximum null distribution.
 */
arma::mat PermTestingIterative(arma::mat data, int nPermutations, int nGroup1);

#endif