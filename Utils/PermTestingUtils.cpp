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
    double colMemMaxMemRatio = totalMemoryPerCol/maxMemoryPerInterval;
    /* Estimate to the nearest hundredth */
    double numPermutationsPerInterval = 100 * floor((1/(100*colMemMaxMemRatio)));
    if((numPermutationsPerInterval*nVoxels*sizeof(double)) > maxMemoryPerInterval)
    {
        std::cout << "Too big of an interval!!!!!!!!!!!" << std::endl;
    }
    return numPermutationsPerInterval;
}