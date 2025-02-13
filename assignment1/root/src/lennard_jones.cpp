#include "lennard_jones.hpp"
#include <cmath>
#include <iostream>


double calculate_sigma_ij(double sigma_i, double sigma_j)
{
    return std::sqrt(sigma_i * sigma_j);
}

double calculate_epsilon_ij(double binding_energy_i, double binding_energy_j)
{
    return std::sqrt(binding_energy_i * binding_energy_j);
}

double calculate_distance(const arma::vec atom_1_coords, const arma::vec atom_2_coords)
{
    arma::vec coords = atom_1_coords - atom_2_coords;

    double dist = std::sqrt(arma::dot(coords, coords));
    return dist;
}

double calculate_lj_energy(double sigma_ij, double radius_ij, double epsilon_ij)
{
    double energy = epsilon_ij * (std::pow((sigma_ij/radius_ij), 12) - 2 * std::pow((sigma_ij/radius_ij),6));
    return energy;
}

double calculate_lj_force(double sigma_ik, double radius_ik, double epsilon_ik)
{
    double analytical_force = epsilon_ik * (1 / radius_ik) * (12 * std::pow(sigma_ik, 12) / std::pow(radius_ik, 13) - 12 * std::pow(sigma_ik, 6) / std::pow(radius_ik, 7));
    return analytical_force;
}
