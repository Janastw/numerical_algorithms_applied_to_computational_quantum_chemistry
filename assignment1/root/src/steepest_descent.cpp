#include "steepest_descent.hpp"
#include <armadillo>

void steepest_descent_central_diff(Cluster cluster, double derivative_step_size, double sd_step_size, double convergence_threshold, int iteration)
{
    if (iteration == 0)
    {
        std::cout << "start steepest descent with golden section line search" << "\n" <<
        "Initial energy: " << cluster.calculate_total_energy() << "\n" <<
        "Stepsize for central difference is:" << derivative_step_size <<
        ";Initial stepsize for line search is:" << sd_step_size <<
        ";Threshold for convergence in force is:" << "" << "\n" <<
        "Central Difference Force" << "\n" << std::endl;
        cluster.update_central_difference(derivative_step_size);
        cluster.print_forward_difference();
    }
    std::cout << "Start steepest descent with golden section line search using central difference force" << "\n" <<
    "Start golden section search\nnew_point" << std::endl;
    cluster.print_forward_difference();
    double original_energy = cluster.calculate_total_energy();
    std::cout << "current energy" << original_energy << std::endl;


    cluster.update_central_difference(derivative_step_size);
    double new_energy = cluster.calculate_total_energy();
    std::cout << "Central Difference Force" << std::endl;
    cluster.print_forward_difference();

    if (new_energy == original_energy)
    {
        std::cout << "Total iterations: " << iteration << "\n" <<
        "Final energy: " << new_energy << "\n" <<
        "Optimized structure:" << std::endl;
        cluster.print_atoms();

        return;
    }
    
    if (new_energy > original_energy)
    {
        steepest_descent_central_diff(cluster, derivative_step_size, sd_step_size / 2, convergence_threshold, iteration + 1);
    }
    else
    {
        steepest_descent_central_diff(cluster, derivative_step_size, sd_step_size * 1.0, convergence_threshold, iteration + 1);
    }
}