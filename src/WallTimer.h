#ifndef _WALL_TIMER_H
#define _WALL_TIMER_H
#include <iostream>
#include <vector>
#include <chrono>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

class WallTimer {
public:
    WallTimer() {
#ifdef __CUDACC__
        cudaEventCreate(&startEvent);
        cudaEventCreate(&stopEvent);
        cudaEventRecord(startEvent, 0);
#else
        start = std::chrono::steady_clock::now();
#endif
        clicks.push_back(start);
    }

    ~WallTimer() {
#ifdef __CUDACC__
        cudaEventDestroy(startEvent);
        cudaEventDestroy(stopEvent);
#endif
    }

    void stop() {
#ifdef __CUDACC__
        cudaEventRecord(stopEvent, 0);
        cudaEventSynchronize(stopEvent);
#else
        end = std::chrono::steady_clock::now();
#endif
    }

    void click() {
#ifdef __CUDACC__
        cudaEvent_t clickEvent;
        cudaEventCreate(&clickEvent);
        cudaEventRecord(clickEvent, 0);
        cudaEventSynchronize(clickEvent);
        clicks.push_back(clickEvent);
#else
        clicks.push_back(std::chrono::steady_clock::now());
#endif
    }

    double elapsed() {
#ifdef __CUDACC__
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, startEvent, stopEvent);
        return milliseconds / 1000.0; // Convert to seconds
#else
        std::chrono::duration<double> elapsed_time = end - start;
        return elapsed_time.count();
#endif
    }

    double elapsedSinceLastClick() {
#ifdef __CUDACC__
        if (clicks.size() < 2) return 0.0f;

        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, clicks[clicks.size() - 2], clicks.back());
        return milliseconds / 1000.0;
#else
        if (clicks.size() < 2) return 0.0;
        std::chrono::duration<double> elapsed_time = clicks.back() - clicks[clicks.size() - 2];
        return elapsed_time.count();
#endif
    }

    double elapsedSinceStart() {
#ifdef __CUDACC__
        if (clicks.empty()) return 0.0f;

        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, startEvent, clicks.back());
        return milliseconds / 1000.0;
#else
        if (clicks.empty()) return 0.0;
        std::chrono::duration<double> elapsed_time = clicks.back() - start;
        return elapsed_time.count();
#endif
    }
    
private:
#ifdef __CUDACC__
    cudaEvent_t startEvent, stopEvent;
    std::vector<cudaEvent_t> clicks;
#else
    std::chrono::steady_clock::time_point start, end;
    std::vector<std::chrono::steady_clock::time_point> clicks;
#endif
};

// Dummy CPU function
void cpu_work(int cycles) {
    for (volatile int i = 0; i < cycles; ++i);
}

#ifdef __CUDACC__
// Dummy CUDA kernel
__global__ void gpu_work(int cycles) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    for (int i = 0; i < cycles; i++);
}
#endif


#endif