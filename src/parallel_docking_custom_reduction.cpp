#include <iostream>
#include <chrono>
#include <vector>
#include <omp.h>

#include "utils.h"

struct BestResult {
    double score;
    int pose_index;
};

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
    auto start_parallel_custom = std::chrono::high_resolution_clock::now();

    BestResult best = {1e100, -1};


    // YOUR CODE HERE


    // End timing
    auto end_parallel_custom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel_custom = end_parallel_custom - start_parallel_custom;

    std::cout << "Interpolated value = " << best.score << std::endl;
    std::cout << "Best pose = " << best.pose_index << std::endl;
    std::cout << "Time taken = " << elapsed_parallel_custom.count() << std::endl;
    return 0;
}
