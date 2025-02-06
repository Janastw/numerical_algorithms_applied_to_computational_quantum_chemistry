#ifndef LENNARD_JONES_HPP
#define LENNARD_JONES_HPP
#include "atom.hpp"
#include "cluster.hpp"


double calculate_sigma_ij(double sigma_i, double sigma_j);

// Binding Energy between i and j
double calculate_epsilon_ij(double binding_energy_i, double binding_energy_j);

// TODO: Can I do const correctness for these functions?
// TODO: Alter this or remove
double calculate_distance(int coord_1, int coord_2);

double calculate_distance(Atom atom_1, Atom atom_2);

double calculate_lennard_jones_energy(Cluster clusters);

void calculate_lennard_jones_forces(Cluster& clusters);

#endif



