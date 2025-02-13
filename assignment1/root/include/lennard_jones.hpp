#ifndef LENNARD_JONES_HPP
#define LENNARD_JONES_HPP
#include "atom.hpp"
#include "cluster.hpp"
#include "derivative_approximation.hpp"
#include <armadillo>



double calculate_sigma_ij(double sigma_i, double sigma_j);

double calculate_epsilon_ij(double binding_energy_i, double binding_energy_j);

double calculate_distance(const arma::vec atom_1_coords, const arma::vec atom_2_coords);

double calculate_lj_energy(double sigma_ij, double radius_ij, double epsilon_ij);

double calculate_lj_force(double sigma_ik, double radius_ik, double epsilon_ik);

#endif



