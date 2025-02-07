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
    // double x = atom_1.x - atom_2.x;
    // double y = atom_1.y - atom_2.y;
    // double z = atom_1.z - atom_2.z;

    arma::vec coords = atom_1.coords - atom_2.coords;

    // double dist = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
    double dist = std::sqrt(arma::dot(coords, coords));
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

    return total_energy;
}


void calculate_lennard_jones_forces(Cluster& clusters)
{
    // TODO: Change clusters.get_num_atoms() to size()?
    int num_atoms = clusters.get_num_atoms();
    std::vector<Atom>& atoms = clusters.get_atoms();

    for (int i = 0; i < num_atoms; i++)
    {
        for (int k = 0; k < num_atoms; k++)
        {
            if (i == k)
            {
                continue;
            }

            double sigma_ik = calculate_sigma_ij(atoms[i].get_sigma(), atoms[k].get_sigma());
            double radius_ik = calculate_distance(atoms[i], atoms[k]);
            double epsilon_ik = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[k].get_epsilon());
            double analytical_force = epsilon_ik * (1 / radius_ik) * (12 * std::pow(sigma_ik, 12) / std::pow(radius_ik, 13) - 12 * std::pow(sigma_ik, 6) / std::pow(radius_ik, 7));
            atoms[i].x_af += analytical_force * (atoms[k].x - atoms[i].x);
            atoms[i].y_af += analytical_force * (atoms[k].y - atoms[i].y);
            atoms[i].z_af += analytical_force * (atoms[k].z - atoms[i].z);
        }
    }
    
}

// void calculate_lennard_jones_forces_forward_difference(Cluster& clusters, double step_size)
// {
//     // Sum of pair-wise energies using sigma_ij/radius_ij
//     int num_atoms = clusters.get_num_atoms();
//     double total_energy = 0;

//     for (int i = 0; i < num_atoms - 1; i++)
//     {
//         for (int j = i + 1; j < num_atoms; j++)
//         {
//             // TODO: Change sigma_ij to be a value retrieved by an atom instance
//             double sigma_ij = calculate_sigma_ij(clusters.get_atoms()[i].get_sigma(), clusters.get_atoms()[j].get_sigma());
//             double radius_ij = calculate_distance(clusters.get_atoms()[i], clusters.get_atoms()[j]);
//             double epsilon_ij = calculate_epsilon_ij(clusters.get_atoms()[i].get_epsilon(), clusters.get_atoms()[j].get_epsilon());
//             total_energy += epsilon_ij * (std::pow((sigma_ij/radius_ij), 12) - 2 * std::pow((sigma_ij/radius_ij),6));
//         }
//     }
    
//     // TODO: Possibly alter the loops to access the atoms with auto instead
//     // for (auto& atom : clusters.get_atoms())
//     // {
//     //     calculate_distance(atom);
//     // }

//     // return total_energy;
// }



/*
Calculate LJ for x, y, and z

for auto& atom: atoms
temp atom.x + step_size

apply step size x, calculate lj
apply step size y, calculate lj
apply step size z, calculate lj

*/
