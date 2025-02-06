#include "cluster.hpp"
#include "lennard_jones.hpp"
#include "derivative_approximation.hpp"
#include <iostream>

int main(void)
{
    // 1. Create cluster
    Cluster gold;
    gold.load_atoms("../sample_input/Force/1.txt");

    // 2. Calculate Energy
    double energy = calculate_lennard_jones_energy(gold);
    gold.print_atoms();

    // 3. Calculate Force
    calculate_lennard_jones_forces(gold);
    gold.print_analytical_force();

    std::cout << "Done" << std::endl;
    return 0;
}
