#include "lennard_jones.hpp"
#include <cmath>
#include <iostream>


double calculate_sigma_ij(double sigma_i, double sigma_j)
{
    return std::sqrt(sigma_i * sigma_j);
}

// Binding Energy between i and j
double calculate_epsilon_ij(double binding_energy_i, double binding_energy_j)
{
    return std::sqrt(binding_energy_i * binding_energy_j);
}

// TODO: Can I do const correctness for these functions?
// TODO: Alter this or remove
double calculate_distance(int coord_1, int coord_2)
{
    return 1;
}

double calculate_distance(Atom atom_1, Atom atom_2)
{
    double x = atom_1.x - atom_2.x;
    double y = atom_1.y - atom_2.y;
    double z = atom_1.z - atom_2.z;

    double dist = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
    return dist;
}

double calculate_lennard_jones_energy(Cluster clusters)
{
    // Sum of pair-wise energies using sigma_ij/radius_ij
    int num_atoms = clusters.get_num_atoms();
    double total_energy = 0;

    for (int i = 0; i < num_atoms - 1; i++)
    {
        for (int j = i + 1; j < num_atoms; j++)
        {
            // TODO: Change sigma_ij to be a value retrieved by an atom instance
            double sigma_ij = calculate_sigma_ij(clusters.get_atoms()[i].get_sigma(), clusters.get_atoms()[j].get_sigma());
            double radius_ij = calculate_distance(clusters.get_atoms()[i], clusters.get_atoms()[j]);
            double epsilon_ij = calculate_epsilon_ij(clusters.get_atoms()[i].get_epsilon(), clusters.get_atoms()[j].get_epsilon());
            total_energy += epsilon_ij * (std::pow((sigma_ij/radius_ij), 12) - 2 * std::pow((sigma_ij/radius_ij),6));
        }
    }
    
    // TODO: Possibly alter the loops to access the atoms with auto instead
    // for (auto& atom : clusters.get_atoms())
    // {
    //     calculate_distance(atom);
    // }

    // TODO: Address chatterbox if necessary
    std::cout << "E_LJ = " << total_energy << std::endl;

    return total_energy;
}


void calculate_lennard_jones_forces(Cluster& clusters)
{
    // TODO: Change clusters.get_num_atoms() to size()?
    int num_atoms = clusters.get_num_atoms();
    std::vector<Atom>& atoms = clusters.get_atoms();

    for (int i = 0; i < num_atoms - 1; i++)
    {
        for (int j = i + 1; j < num_atoms; j++)
        {
            // TODO: Change sigma_ij to be a value retrieved by an atom instance
            double sigma_ij = calculate_sigma_ij(atoms[i].get_sigma(), atoms[j].get_sigma());
            double radius_ij = calculate_distance(atoms[i], atoms[j]);
            double epsilon_ij = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[j].get_epsilon());
            // double energy = calc.calculate_energy(clusters);
            double new_force = 12 * epsilon_ij / (radius_ij) * (std::pow((sigma_ij/radius_ij), 12) - std::pow((sigma_ij/radius_ij),6));
            atoms[i].x_af += new_force * (atoms[j].x - atoms[i].x);
            atoms[i].y_af += new_force * (atoms[j].y - atoms[i].y);
            atoms[i].z_af += new_force * (atoms[j].z - atoms[i].z);
        }
    }
}
