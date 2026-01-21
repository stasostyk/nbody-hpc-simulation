#pragma once

#include <chrono>
#include <iostream>

class Timer {
public:
    Timer() {}
    ~Timer() = default;

    void start() {
        tStart = std::chrono::steady_clock::now();
    }

    void end() {
        tEnd = std::chrono::steady_clock::now();
    }

    void print() {
        std::cout << "Ellapsed time: " << ellapsedMs() << " ms" << std::endl;
    }

    double ellapsedMs() {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart);
        return (double) duration.count();
    }

    std::chrono::steady_clock::time_point tStart;
    std::chrono::steady_clock::time_point tEnd;

};