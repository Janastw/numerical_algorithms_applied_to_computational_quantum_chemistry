#include "cluster.hpp"
#include "lennard_jones.hpp"
#include "derivative_approximation.hpp"
#include <iostream>
#include <fstream>

void process_energy_file(int file_number) {
    std::string input_filename = "./sample_input/Energy/" + std::to_string(file_number) + ".txt";
    std::string output_filename = "project_outputs/Energy_" + std::to_string(file_number) + ".txt";
    std::filesystem::create_directories("project_outputs");

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Error: Could not open " << output_filename << " for writing.\n";
        return;
    }

    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outfile.rdbuf());

    Cluster cluster;
    if (!cluster.load_atoms(input_filename)) {
        std::cerr << "Error: Failed to load " << input_filename << std::endl;
        return;
    }

    double energy = cluster.calculate_total_energy();
    cluster.print_atoms();
    std::cout << "E_LJ = " << energy << "\n" << std::endl;

    std::cout.rdbuf(coutbuf);
    outfile.close();
}

void process_standard_SD_file(int file_number) {
    std::string input_filename = "./sample_input/standard_SD/" + std::to_string(file_number) + ".txt";
    std::string output_filename = "project_outputs/standard_SD_" + std::to_string(file_number) + ".txt";

    std::filesystem::create_directories("project_outputs");

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Error: Could not open " << output_filename << " for writing.\n";
        return;
    }

    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outfile.rdbuf());

    Cluster gold_SD;
    if (!gold_SD.load_atoms(input_filename)) {
        std::cerr << "Error: Failed to load " << input_filename << std::endl;
        return;
    }

    double energy_standard_SD = gold_SD.calculate_total_energy();
    double deriv_approx_step_size = 0.0001;
    double sd_step_size = 0.3;
    double threshold_convergence = 0.01;

    std::cout << "Start steepest descent with golden section line search\n"
              << "Initial energy: " << energy_standard_SD << "\n"
              << "Stepsize for central difference is: " << deriv_approx_step_size << "\n"
              << "Initial stepsize for line search is: " << sd_step_size << "\n"
              << "Threshold for convergence in force is: " << threshold_convergence << "\n"
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
    
    std::cout << "Final energy: " << gold_SD.calculate_total_energy() << std::endl;
    std::cout << "Optimized structure:\n" << gold_SD.system_coordinates << std::endl;
    
    std::cout.rdbuf(coutbuf);
    outfile.close();
}

void process_force_file(int file_number) {
    std::string input_filename = "./sample_input/Force/" + std::to_string(file_number) + ".txt";
    std::string output_filename = "project_outputs/Force_" + std::to_string(file_number) + ".txt";

    std::filesystem::create_directories("project_outputs");

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Error opening " << output_filename << " for writing." << std::endl;
        return;
    }

    std::streambuf* coutbuf = std::cout.rdbuf();  // Save original buffer
    std::cout.rdbuf(outfile.rdbuf());  // Redirect std::cout to file

    Cluster gold_force;
    if (!gold_force.load_atoms(input_filename)) {
        std::cerr << "Error: Failed to load " << input_filename << std::endl;
        return;
    }

    double energy_force = gold_force.calculate_total_energy();
    std::cout << "E_LJ = " << energy_force << std::endl;
    
    gold_force.update_analytical_force();
    gold_force.print_analytical_force();
    // gold_force.print_system_coordinate_forces();
    gold_force.zero_forces();

    double step_sizes[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001};
    for (double step : step_sizes) {
        std::cout << "Stepsize for finite difference: " << step << std::endl;
        
        std::cout << "F_LJ_forward_difference" << std::endl;
        gold_force.update_forward_difference(step);
        gold_force.print_analytical_force();
        // gold_force.print_system_coordinate_forces();
        gold_force.zero_forces();

        std::cout << "F_LJ_central_difference" << std::endl;
        gold_force.update_central_difference(step);
        gold_force.print_analytical_force();
        // gold_force.print_system_coordinate_forces();
        gold_force.zero_forces();
    }

    std::cout.flush();
    std::cout.rdbuf(coutbuf);
    outfile.close();
}

void process_SD_with_line_search(int file_number) {
    std::string input_filename = "./sample_input/SD_with_line_search/" + std::to_string(file_number) + ".txt";
    std::string output_filename = "project_outputs/SD_with_line_search_" + std::to_string(file_number) + ".txt";

    std::filesystem::create_directories("project_outputs");

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Error opening " << output_filename << " for writing." << std::endl;
        return;
    }

    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outfile.rdbuf());

    Cluster gold_SD_with_line_search;
    if (!gold_SD_with_line_search.load_atoms(input_filename)) {
        std::cerr << "Error: Failed to load " << input_filename << std::endl;
        return;
    }

    double energy_standard_SD = gold_SD_with_line_search.calculate_total_energy();
    double deriv_approx_step_size = 0.0001;
    double sd_step_size = 0.3;
    double threshold_convergence = 0.01;

    std::cout << "Start steepest descent with golden section line search\n"
              << "Initial energy: " << energy_standard_SD
              << "\nStepsize for central difference is: " << deriv_approx_step_size
              << "\nInitial stepsize for line search is: " << sd_step_size
              << "\nThreshold for convergence in force is: " << threshold_convergence
              << "\nAnalytical Force" << std::endl;

    gold_SD_with_line_search.update_analytical_force();
    gold_SD_with_line_search.print_system_coordinate_forces();
    
    std::cout << "Forward Difference Force" << std::endl;
    gold_SD_with_line_search.zero_forces();
    gold_SD_with_line_search.update_forward_difference(deriv_approx_step_size);
    gold_SD_with_line_search.print_system_coordinate_forces();
    
    std::cout << "Central Difference Force" << std::endl;
    gold_SD_with_line_search.zero_forces();
    gold_SD_with_line_search.update_central_difference(deriv_approx_step_size);
    gold_SD_with_line_search.print_system_coordinate_forces();

    gold_SD_with_line_search.steepest_descent(deriv_approx_step_size, sd_step_size, threshold_convergence);
    
    std::cout << "Final energy: " << gold_SD_with_line_search.calculate_total_energy() << std::endl;
    std::cout << gold_SD_with_line_search.system_coordinates << std::endl;

    std::cout.flush();
    std::cout.rdbuf(coutbuf);
    outfile.close();
}

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

    // 1. Calculate Energy
    // If files do not output, uncomment the above section to print in terminal
    std::vector<int> energy_files = {1, 2, 3};

    for (int file_number : energy_files) {
        process_energy_file(file_number);
    }

    // 2. Calculate Force

    // Cluster gold_force_1;
    // gold_force_1.load_atoms("./sample_input/Force/1.txt");
    // double energy_force_1 = gold_force_1.calculate_total_energy();
    // std::cout << "E_LJ = " << energy_force_1 << std::endl;
    // gold_force_1.update_analytical_force();
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.1);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.1);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.01);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.01);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // std::cout << "Stepsize for finite difference:" << std::endl;
    // std::cout << "F_LJ_forward_difference" << std::endl;
    // gold_force_1.update_forward_difference(0.001);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();
    // std::cout << "F_LJ_central_difference" << std::endl;
    // gold_force_1.update_central_difference(0.001);
    // gold_force_1.print_analytical_force();
    // // gold_force_1.print_system_coordinate_forces();
    // gold_force_1.zero_forces();

    // Question 2 - Calculate Force
    // If File does not output, please uncomment the code above and it will print in the terminal
    std::vector<int> force_file_numbers = {1};  // Add more numbers if needed

    for (int file_number : force_file_numbers) {
        process_force_file(file_number);
    }

    // Question 3 - standard_SD

    // Cluster gold_SD;
    // gold_SD.load_atoms("./sample_input/standard_SD/1.txt");
    // double energy_standard_SD = gold_SD.calculate_total_energy();
    // double deriv_approx_step_size = 0.0001;
    // double sd_step_size = 0.3;
    // double threshold_convergence = 0.01;

    // std::cout << "Start steepest descent with golden section line search\n"
    // << "Initial energy: " << energy_standard_SD
    // << "\nStepsize for central difference is:" << deriv_approx_step_size
    // << "Initial stepsize for line search is:" << sd_step_size
    // << "Threshold for convergence in force is:" << threshold_convergence
    // << "Analytical Force" << std::endl;
    // gold_SD.update_analytical_force();
    // gold_SD.print_system_coordinate_forces();
    // std::cout << "Forward Difference Force" << std::endl;
    // gold_SD.zero_forces();
    // gold_SD.update_forward_difference(deriv_approx_step_size);
    // gold_SD.print_system_coordinate_forces();
    // std::cout << "Central Difference Force" << std::endl;
    // gold_SD.zero_forces();
    // gold_SD.update_central_difference(deriv_approx_step_size);
    // gold_SD.print_system_coordinate_forces();

    // gold_SD.steepest_descent(deriv_approx_step_size, sd_step_size, threshold_convergence);
    // std::cout << "Final energy:" << gold_SD.calculate_total_energy() << std::endl;
    // std::cout << gold_SD.system_coordinates << std::endl;

    std::vector<int> file_numbers = {1};

    for (int file_number : file_numbers) {
        process_standard_SD_file(file_number);
    }

    // Question 3 - SD with line search

    // Cluster gold_SD_with_line_search;
    // gold_SD_with_line_search.load_atoms("./sample_input/SD_with_line_search/1.txt");
    // double energy_standard_SD = gold_SD_with_line_search.calculate_total_energy();
    // double deriv_approx_step_size = 0.0001;
    // double sd_step_size = 0.3;
    // double threshold_convergence = 0.01;

    // std::cout << "Start steepest descent with golden section line search\n"
    // << "Initial energy: " << energy_standard_SD
    // << "\nStepsize for central difference is:" << deriv_approx_step_size
    // << "Initial stepsize for line search is:" << sd_step_size
    // << "Threshold for convergence in force is:" << threshold_convergence
    // << "Analytical Force" << std::endl;
    // gold_SD_with_line_search.update_analytical_force();
    // gold_SD_with_line_search.print_system_coordinate_forces();
    // std::cout << "Forward Difference Force" << std::endl;
    // gold_SD_with_line_search.zero_forces();
    // gold_SD_with_line_search.update_forward_difference(deriv_approx_step_size);
    // gold_SD_with_line_search.print_system_coordinate_forces();
    // std::cout << "Central Difference Force" << std::endl;
    // gold_SD_with_line_search.zero_forces();
    // gold_SD_with_line_search.update_central_difference(deriv_approx_step_size);
    // gold_SD_with_line_search.print_system_coordinate_forces();

    // gold_SD_with_line_search.steepest_descent(deriv_approx_step_size, sd_step_size, threshold_convergence);
    // std::cout << "Final energy:" << gold_SD_with_line_search.calculate_total_energy() << std::endl;
    // std::cout << gold_SD_with_line_search.system_coordinates << std::endl;

    // Question 3 - SD with line search
    // If files don't generate, please uncomment above and view the outputs in the terminal
    
    std::vector<int> line_file_numbers = {1, 2};

    for (int file_number : line_file_numbers) {
        process_SD_with_line_search(file_number);
    }

    return 0;
}
