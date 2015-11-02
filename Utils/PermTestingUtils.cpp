#include <math.h>
#include <iostream>
/**
 * @brief GetIntervalDimension
 * @details 
 *         See PermTesting.h for information
 */
int GetIntervalDimension(int nVoxels, double maxMemoryPerInterval)
{
    double totalMemoryPerCol = (double)(nVoxels * sizeof(double));
    double numPermutationsPerInterval = maxMemoryPerInterval/totalMemoryPerCol;
    /* Estimate to the nearest hundredth */
    int intervalLength = floor(numPermutationsPerInterval);

    return intervalLength;
}