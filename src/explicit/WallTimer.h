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
        clicks.push_back(startEvent);  // Use correct type
#else
        start = std::chrono::steady_clock::now();
        clicks.push_back(start);
#endif
    }

    ~WallTimer() {
#ifdef __CUDACC__
        cudaEventDestroy(startEvent);
        cudaEventDestroy(stopEvent);
        for (auto& e : clicks) cudaEventDestroy(e);
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

    double elapsed() const {
#ifdef __CUDACC__
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, startEvent, stopEvent);
        return milliseconds / 1000.0;
#else
        return std::chrono::duration<double>(end - start).count();
#endif
    }

    double elapsedSinceLastClick() const {
        if (clicks.size() < 2) return 0.0;
#ifdef __CUDACC__
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, clicks[clicks.size() - 2], clicks.back());
        return milliseconds / 1000.0;
#else
        return std::chrono::duration<double>(
            clicks.back() - clicks[clicks.size() - 2]
        ).count();
#endif
    }

    double elapsedSinceStart() const {
        if (clicks.empty()) return 0.0;
#ifdef __CUDACC__
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, startEvent, clicks.back());
        return milliseconds / 1000.0;
#else
        return std::chrono::duration<double>(
            clicks.back() - start
        ).count();
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
inline void cpu_work(int cycles) {
    for (volatile int i = 0; i < cycles; ++i);
}

#ifdef __CUDACC__
// Dummy CUDA kernel
__global__ void gpu_work(int cycles) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    for (int i = 0; i < cycles; i++);
}
#endif

#endif // _WALL_TIMER_H
