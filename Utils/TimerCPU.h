#ifndef TIMER_CPU_H
#define TIMER_CPU_H

// --------------------------------------------------------------------
// CpuTimer
//
// This utility class encapsulates a simple timer for recording the
// time between a start and stop event.
// --------------------------------------------------------------------
//  Windows
#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#endif

#include <string>
#include <iostream>

class CpuTimer {
public:
    CpuTimer() : m_start(0), m_stop(0) {}
    CpuTimer(std::string TimerName) : m_start(0), m_stop(0) 
    {
        name = TimerName;
    }
    ~CpuTimer() {}

    // wall time
    void Start() { m_start = get_wall_time(); }
    void Stop()  { m_stop = get_wall_time(); }
    float Elapsed() {
        return (m_stop - m_start);
    }

    // cpu time
    void Start_cputimer() { m_start_cpu = get_cpu_time(); }
    void Stop_cputimer()  { m_stop_cpu = get_cpu_time(); }
    float Elapsed_cputimer() {
        return (m_stop_cpu - m_start_cpu);
    }

    void print()
    {
        std::cout << name <<  " wall time = " << Elapsed() << " seconds" << std::endl;
        std::cout << name <<  " cpu time = " << Elapsed_cputimer() << " seconds" << std::endl;
    }

private:
#ifdef _WIN32
    float get_wall_time(){
        LARGE_INTEGER time, freq;
        if (!QueryPerformanceFrequency(&freq)){
            //  Handle error
            return 0;
        }
        if (!QueryPerformanceCounter(&time)){
            //  Handle error
            return 0;
        }
        return (float)time.QuadPart / freq.QuadPart;
    }

    float get_cpu_time(){
        FILETIME a, b, c, d;
        if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0){
            //  Returns total user time.
            //  Can be tweaked to include kernel times as well.
            return
                (float)(d.dwLowDateTime |
                ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
        }
        else{
            //  Handle error
            return 0;
        }
    }
#else
    #if OPENMP_ENABLED
    float get_wall_time()
    {
        return omp_get_wtime();
    }
    #else
    float get_wall_time()
    {
        struct timeval time;
        if (gettimeofday(&time, NULL)){
            //  Handle error
            return 0;
        }
        return (float)time.tv_sec + (float)time.tv_usec * .000001;
    }
    #endif


    float get_cpu_time(){
        return (double)clock() / CLOCKS_PER_SEC;
    }
#endif
private:
    // wall time
    float m_start;
    float m_stop;

    // cpu time
    float m_start_cpu;
    float m_stop_cpu;

    std::string name;
};


#endif
