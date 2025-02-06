#include "cluster.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

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

void Cluster::print_analytical_force()
{
    if (atoms.empty())
    {
        std::cout << "No atoms present in the system" << std::endl;
        return;
    }
    for (int i = 0; i < atoms.size(); i ++)
    {
        std::cout << atoms[i].x_af << " ";
    }
        std::cout << std::endl;
    for (int i = 0; i < atoms.size(); i ++)
    {
        std::cout << atoms[i].y_af << " ";
    }
        std::cout << std::endl;
    for (int i = 0; i < atoms.size(); i ++)
    {
        std::cout << atoms[i].z_af << " ";
    }
    std::cout << std::endl;
    for (const auto& atom : atoms)
    {
        // std::cout << atom.x << " " << atom.y << " " << atom.z << std::endl;
        // std::cout << atom.x_af << " " << atom.y_af << " " << atom.z_af << std::endl;
    }
}