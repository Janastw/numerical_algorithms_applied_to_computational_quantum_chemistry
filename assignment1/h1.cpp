// #include <armadillo>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>
#include <cmath>


// TODO: Raise an error if the atomic number is not 79 (Gold/Au) - use std::exception
// Echo the input

class Atom
{
    private:
        int atomic_number;
        std::string bond_units = "angstroms";
        double sigma;
    
    public:
        // TODO Use eigen or armadillo for coordinates, but we will worry about that later
        double x;
        double y;
        double z;

        Atom(int atomic_number_, int x_, int y_, int z_) :
            atomic_number(atomic_number_), x(x_), y(y_), z(z_)
            {
                if (atomic_number_ != 79) // Error handles non-gold. TODO handle < 1 cases
                {
                    throw std::exception();
                }
                sigma = 0.26290;
            }
        
        int get_atomic_number() const { atomic_number; }

        // TODO Change to armadillo vector
        void get_coords() const {}

        double get_sigma() const { sigma; }

        void set_atomic_number(int new_atomic_number)
        {
            atomic_number = new_atomic_number;
        }

        void set_coords()
        {

        }

        Atom operator-(Atom atom_2)
        {
            x - atom_2.x;
        }
};

// TODO Getter and Setter Functions
class Cluster
{
    private:
        std::vector<Atom> atoms;
        int num_atoms = 0;

    public:

        // TODO: Find out if I have to template this or if I can just take in an integer, float, double etc.
        void add_atom(int atomic_number, double x, double y, double z)
        {
            Atom new_atom(atomic_number, x, y, z);
            atoms.push_back(new_atom);
            num_atoms++;   
        }

        int get_num_atoms() const { return num_atoms; }

        std::vector<Atom> get_atoms() const { return atoms; }

        void print_atoms()
        {
            std::cout << num_atoms << std::endl;

            for (auto& atom : atoms)
            {
                std::cout << atom.get_atomic_number() << "(" << ")" << std::endl;
            }
        }
};

// Load atoms in a system
// Check boundaries?


// TODO Change to allow file to be read through command line
    // TODO Allow for an entire folder to be run too
Cluster load_atoms(std::string file)
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

    // Echoing input to output
    std::cout << line << std::endl;;

    Cluster cluster_1;

    // Input files have atomic number, x, y, then z coordinate in each line
    while (getline(inputFile, line))
    {
        // Echoing input to output
        std::cout << line << std::endl;

        int atomic_number;
        double x;
        double y;
        double z;
        
        std::istringstream numbers(line);
        numbers >> atomic_number;
        numbers >> x;
        numbers >> y;
        numbers >> z;

        cluster_1.add_atom(atomic_number, x, y, z);
        
    }

    return cluster_1;
}

// TODO: Can I do const correctness for these functions?
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
    std::cout << "printing distance: " << dist << std::endl;
    return dist;
}

double calculate_energy(Cluster clusters)
{
    // Sum of pair-wise energies using sigma_ij/radius_ij
    int num_atoms = clusters.get_num_atoms();
    double total_energy = 0;

    for (int i = 0; i < num_atoms - 1; i++)
    {
        for (int j = i + 1; j < num_atoms; j++)
        {
            double sigma_ij = 0.2629;
            double radius_ij = calculate_distance(clusters.get_atoms()[i], clusters.get_atoms()[j]);
            total_energy += std::pow((sigma_ij/radius_ij), 6) - 2 * std::pow((sigma_ij/radius_ij), 12);
            std::cout << total_energy << std::endl;
        }
    }
    
    // for (auto& atom : clusters.get_atoms())
    // {
    //     calculate_distance(atom);
    // }
    
    return total_energy;
}

int main(void)
{
    Cluster gold = load_atoms("./sample_input/Energy/1.txt");
    // gold.print_atoms();

    // 1. Calculate Energy
    double energy_of_the_cluster = calculate_energy(gold); // calculate_energy();
    std::cout << "E_LJ = " << energy_of_the_cluster << std::endl;

    // 2. Calculate Force
    int force_of_the_cluster = 0;
    std::cout << "" << std::endl;


    std::cout << "Done" << std::endl;
    return 0;
}