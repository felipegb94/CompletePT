#ifndef PERMTESTINGUTILS_H
#define PERMTESTINGUTILS_H
/**
 * @brief GetIntervalDimension
 * @details 
 *         Calculates the number of permutations that can be done given the maximum amount of 
 *         memory that can be used at a time and the length of the output of a permutation 
 *         (nVoxel).
 * @param nVoxels Number of voxels per subject. This is the output of doing one permutation.
 * @param maxMemoryPerInterval. Maximum memory used per matrix multiplication. Value given in bytes
 * . For example 5e8 (500 megabytes), 1e9 (1 gigabyte).
 * @return numPermutationsPerInterval  - How many permutation can be done given the 
 * maxMemoryPerInterval
 */
int GetIntervalDimension(int nVoxels, double maxMemoryPerInterval);

#endif