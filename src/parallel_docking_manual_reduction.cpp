#include <iostream>
#include <chrono>
#include <vector>
#include <omp.h>

#include "utils.h"

int main()
{
    Grid grid = read_grid("../data/grid.pts");
    std::vector<Atom> ligand_atoms = read_xyz("../data/ligand.xyz");

    // Generate poses (serial, unchanged)
    std::vector<std::vector<Atom>> poses;
    poses.reserve(1000000);
    for (int i = 0; i < 1000000; i++) {
        poses.push_back(transform_ligand(ligand_atoms, i));
    }

    // Start timing
    auto start_parallel = std::chrono::high_resolution_clock::now();

    double global_min = 10000;
    int best_pose = -1;


    // YOUR CODE HERE


    // End timing
    auto end_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel = end_parallel - start_parallel;

    std::cout << "Interpolated value = " << global_min << std::endl;
    std::cout << "Best pose = " << best_pose << std::endl;
    std::cout << "Time taken = " << elapsed_parallel.count() << std::endl;

    return 0;
}
