#include "cluster.hpp"
#include "lennard_jones.hpp"
#include "derivative_approximation.hpp"
#include "steepest_descent.hpp"
#include <iostream>
#include <fstream>

int main(void)
{
    // 1. Calculate Energy
    // Cluster gold_energy_1;
    // gold_energy_1.load_atoms("./sample_input/Energy/1.txt");
    // double energy_1 = gold_energy_1.calculate_total_energy();
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
    // Cluster gold_force_1;
    // gold_force_1.load_atoms("./sample_input/Force/1.txt");
    // double energy_force_1 = gold_force_1.calculate_total_energy();
    // std::cout << "E_LJ = " << energy_force_1 << std::endl;
    // gold_force_1.update_analytical_force();
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.1);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.1);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.01);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.01);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.001);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.001);
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // FILE OUTPUT FOR QUESTION 2 PORTION

    // std::ofstream outfile("project_outputs/output.txt");
    // if (!outfile) {
    //     std::cerr << "Error opening file for writing." << std::endl;
    //     return 1;
    // }

    // std::streambuf* coutbuf = std::cout.rdbuf(); // Save old buffer
    // std::cout.rdbuf(outfile.rdbuf()); // Redirect std::cout to file

    // Cluster gold_force_1;
    // gold_force_1.load_atoms("./sample_input/Force/1.txt");
    // double energy_force_1 = gold_force_1.calculate_total_energy();
    // std::cout << "E_LJ = " << energy_force_1 << std::endl;
    // gold_force_1.update_analytical_force();
    // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // double step_sizes[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001};
    // for (double step : step_sizes) {
    //     std::cout << "Stepsize for finite difference: " << step << std::endl;
        
    //     std::cout << "F_LJ_forward_difference" << std::endl;
    //     gold_force_1.update_forward_difference(step);
    //     gold_force_1.print_system_coordinate_forces();
    //     gold_force_1.zero_forces();

    //     std::cout << "F_LJ_central_difference" << std::endl;
    //     gold_force_1.update_central_difference(step);
    //     gold_force_1.print_system_coordinate_forces();
    //     gold_force_1.zero_forces();
    // }

    // std::cout.flush(); // Ensure everything is written
    // std::cout.rdbuf(coutbuf); // Restore old buffer
    // outfile.close(); // Close file properly

    // std::cout << "Output successfully written to output.txt" << std::endl;


    // Question 3 - standard_SD

    Cluster gold_SD;
    gold_SD.load_atoms("./sample_input/standard_SD/1.txt");
    double energy_standard_SD = gold_SD.calculate_total_energy(gold_SD.system_coordinates);
    double deriv_approx_step_size = 0.0001;
    double sd_step_size = 0.3;
    double threshold_convergence = 0.01;

    std::cout << "Start steepest descent with golden section line search\n"
    << "Initial energy: " << energy_standard_SD
    << "\nStepsize for central difference is:" << deriv_approx_step_size
    << "Initial stepsize for line search is:" << sd_step_size
    << "Threshold for convergence in force is:" << threshold_convergence
    << "Analytical Force" << std::endl;
    gold_SD.update_analytical_force();
    gold_SD.print_system_coordinate_forces();
    std::cout << "Forward Difference Force" << std::endl;
    gold_SD.zero_forces();
    gold_SD.update_forward_difference(deriv_approx_step_size);
    gold_SD.print_system_coordinate_forces();
    std::cout << "Central Difference Force" << std::endl;
    gold_SD.zero_forces();
    gold_SD.update_central_difference(deriv_approx_step_size);
    gold_SD.print_system_coordinate_forces();

    gold_SD.steepest_descent(deriv_approx_step_size, sd_step_size, threshold_convergence);
    std::cout << "Final energy:" << gold_SD.calculate_total_energy(gold_SD.system_coordinates) << std::endl;
    std::cout << gold_SD.system_coordinates << std::endl;


    // gold_SD.update_analytical_force();
    // std::cout << gold_SD.system_coordinate_forces << std::endl;
    // gold_SD.zero_forces();
    // gold_SD.update_forward_difference(0.0001);
    // std::cout << gold_SD.system_coordinate_forces << std::endl;
    // gold_SD.zero_forces();
    // gold_SD.update_central_difference(0.0001);
    // std::cout << gold_SD.system_coordinate_forces << std::endl;
    // gold_SD.zero_forces();
    // std::cout << "Final energy: " << gold_SD.calculate_total_energy() << std::endl;
    // std::cout << "Optimized structure: " << std::endl;
    // std::cout << gold_SD.system_coordinates << std::endl;



    std::cout << "\n" <<"Done" << std::endl;
    return 0;
}
