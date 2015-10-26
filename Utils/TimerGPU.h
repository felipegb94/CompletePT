#ifndef TIMER_GPU_H
#define TIMER_GPU_H

#include <string>
#include <iostream>
// ----------------------------------------------------------------------------
// CUDA headers
// ----------------------------------------------------------------------------
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_launch_parameters.h>


// --------------------------------------------------------------------
// GpuTimer
//
// This utility class encapsulates a simple timer for recording the
// time between a start and stop event.
// --------------------------------------------------------------------
class GpuTimer
{
public:
    GpuTimer(cudaStream_t stream = 0) : m_stream(stream) {
        cudaEventCreate(&m_start);
        cudaEventCreate(&m_stop);
    }
    GpuTimer(std::string timerName, cudaStream_t stream = 0) : m_stream(stream) {
        cudaEventCreate(&m_start);
        cudaEventCreate(&m_stop);
        name = timerName;
    }

    ~GpuTimer() {
        cudaEventDestroy(m_start);
        cudaEventDestroy(m_stop);
    }

    void Start() {cudaEventRecord(m_start, m_stream);}
    void Stop()  {cudaEventRecord(m_stop,  m_stream);}
    void print()
    {
        std::cout << name <<  " execution time = " << Elapsed()/1000 << " seconds" << std::endl;
    }

    float Elapsed() {
        float elapsed;
        cudaEventSynchronize(m_stop);
        cudaEventElapsedTime(&elapsed, m_start, m_stop);
        return elapsed;
    }


private:
    std::string name;
    cudaStream_t m_stream;
    cudaEvent_t  m_start;
    cudaEvent_t  m_stop;
};

#endif
