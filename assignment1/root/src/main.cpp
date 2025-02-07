#include "cluster.hpp"
#include "lennard_jones.hpp"
#include "derivative_approximation.hpp"
#include <iostream>

int main(void)
{
    // 1. Calculate Energy
    // Cluster gold_energy_1;
    // gold_energy_1.load_atoms("./sample_input/Energy/1.txt");
    // double energy_1 = gold_energy_1.calculate_total_energy();
    // // double energy_1 = calculate_lennard_jones_energy(gold_energy_1);
    // gold_energy_1.print_atoms();
    // std::cout << "E_LJ = " << energy_1 << "\n" << std::endl;

    // Cluster gold_energy_2;
    // gold_energy_2.load_atoms("./sample_input/Energy/2.txt");
    // double energy_2 = gold_energy_2.calculate_total_energy();
    // gold_energy_2.print_atoms();
    // std::cout << "E_LJ = " << energy_2 << "\n" << std::endl;

    // Cluster gold_energy_3;
    // gold_energy_3.load_atoms("./sample_input/Energy/3.txt");
    // double energy_3 = gold_energy_3.calculate_total_energy();
    // gold_energy_3.print_atoms();
    // std::cout << "E_LJ = " << energy_3 << std::endl;

    // 2. Calculate Force
    Cluster gold_force_1;
    gold_force_1.load_atoms("./sample_input/Force/1.txt");
    double energy_force_1 = gold_force_1.calculate_total_energy();
    std::cout << "E_LJ = " << energy_force_1 << std::endl;
    // calculate_lennard_jones_forces(gold_force_1);
    gold_force_1.update_analytical_force();
    gold_force_1.print_analytical_force();
    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // std::cout << calculate_lennard_jones_forces_forward_difference(gold_force_1, 0.1) << std::endl;

    std::cout << "\n" <<"Done" << std::endl;
    return 0;
}
