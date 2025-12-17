#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>

// -----------------------------------------
// Serial vector norm
// -----------------------------------------
double norm_serial(const std::vector<double>& v) {
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return std::sqrt(sum);
}

// -----------------------------------------
// Parallel vector norm
// -----------------------------------------
double norm_parallel(const std::vector<double>& v) {
    double sum = 0.0;

    // YOUR CODE HERE

    return std::sqrt(sum);
}

int main() {
    // Make vector size ~50,000,000 elements (~400 MB of memory)
    const size_t N = 50000000;
    std::vector<double> v(N, 1.0);  // Fill vector with 1.0's

    // Serial timing
    auto start_serial = std::chrono::high_resolution_clock::now();
    double out_serial = norm_serial(v);
    auto end_serial = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_serial = end_serial - start_serial;
    double t_serial = elapsed_serial.count();

    // Parallel timing
    auto start_parallel = std::chrono::high_resolution_clock::now();
    double out_parallel = norm_parallel(v);
    auto end_parallel = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_parallel = end_parallel - start_parallel;
    double t_parallel = elapsed_parallel.count();

    // Print
    std::cout << "Serial norm    = " << out_serial << "\n";
    std::cout << "Parallel norm  = " << out_parallel << "\n\n";

    std::cout << "Serial time    = " << t_serial    << " s\n";
    std::cout << "Parallel time  = " << t_parallel  << " s\n";
    std::cout << "Speedup        = " << t_serial / t_parallel << "x\n";

    return 0;
}