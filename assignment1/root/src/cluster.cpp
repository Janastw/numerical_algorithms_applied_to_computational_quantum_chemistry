#include "cluster.hpp"
#include <fstream>
#include <iostream>
#include <sstream>


// Cluster::Cluster(const Cluster& cluster_, double step_size) :
//     atoms(cluster_.get_atoms())
//     {
//         for (auto& atom : atoms)
//         {
//             atom.x += step_size;
//             atom.y += step_size;
//             atom.z += step_size;
//         }
//     }

bool Cluster::load_atoms(std::string file)
{
    std::ifstream inputFile(file);
    if (!inputFile)
    {
        // TODO: Output error message
        // Incorrect or corrupt file
        throw std::exception();
    }
    std::string line;
    getline(inputFile, line);

    // First line is the number of atoms
    int num_atoms = std::stoi(line);

    // Input files have atomic number, x, y, then z coordinate in each line
    while (getline(inputFile, line))
    {

        int atomic_number;
        double x;
        double y;
        double z;
        
        std::istringstream numbers(line);
        numbers >> atomic_number;
        numbers >> x;
        numbers >> y;
        numbers >> z;

        add_atom(atomic_number, x, y, z);
        
    }
    return true;
}

// TODO: Find out if I have to template this or if I can just take in an integer, float, double etc.
void Cluster::add_atom(int atomic_number, double x, double y, double z)
{
    Atom new_atom(atomic_number, x, y, z);
    atoms.push_back(new_atom);
}

int Cluster::get_num_atoms() const { return atoms.size(); }

std::vector<Atom> Cluster::get_atoms() const { return atoms; }
std::vector<Atom>& Cluster::get_atoms() { return atoms; }

void Cluster::print_atoms()
{
    for (auto& atom : atoms)
    {
        std::cout << atom.get_atomic_number() << "(" << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }
}

double Cluster::calculate_total_energy()
{
    // Sum of pair-wise energies using sigma_ij/radius_ij
    int num_atoms = atoms.size();
    double total_energy = 0;

    for (int i = 0; i < num_atoms - 1; i++)
    {
        for (int j = i + 1; j < num_atoms; j++)
        {
            // TODO: Change sigma_ij to be a value retrieved by an atom instance
            double sigma_ij = calculate_sigma_ij(atoms[i].get_sigma(), atoms[j].get_sigma());
            double radius_ij = calculate_distance(atoms[i], atoms[j]);
            double epsilon_ij = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[j].get_epsilon());
            total_energy += calculate_lj_energy(sigma_ij, radius_ij, epsilon_ij);
        }
    }
    
    // TODO: Possibly alter the loops to access the atoms with auto instead
    // for (auto& atom : clusters.get_atoms())
    // {
    //     calculate_distance(atom);
    // }

    return total_energy;
}

void Cluster::update_analytical_force()
{
    int num_atoms = atoms.size();
    std::vector<Atom>& atoms = get_atoms();

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
            double analytical_force = calculate_lj_force(sigma_ik, radius_ik, epsilon_ik);

            atoms[i].coords_analytical_forces += (atoms[i].coords - atoms[k].coords) * analytical_force;
        }
    }
}

void Cluster::print_analytical_force()
{
    if (atoms.empty())
    {
        std::cout << "No atoms present in the system" << std::endl;
        return;
    }
    for (int i = 0; i < 3; i++)
    {
        for (auto& atom: atoms)
        {
            std::cout << std::setw(8) << std::setprecision(4) << atom.coords_analytical_forces[i] << " ";
        }
        std::cout << std::endl;
    }

}

void Cluster::print_forward_difference()
{
    if (atoms.empty())
    {
        std::cout << "No atoms present in the system" << std::endl;
        return;
    }
    if (atoms.empty())
    {
        std::cout << "No atoms present in the system" << std::endl;
        return;
    }
    for (int i = 0; i < 3; i++)
    {
        for (auto& atom: atoms)
        {
            std::cout << std::setw(8) << std::setprecision(4) << atom.coords_analytical_forces[i] << " ";
        }
        std::cout << std::endl;
    }
}