// #include <armadillo>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>
#include <cmath>
#include <armadillo>


// TODO: Raise an error if the atomic number is not 79 (Gold/Au) - use std::exception
// Echo the input

class Atom
{
    private:
        int atomic_number;
        std::string bond_units = "angstroms";
        double sigma; // angstroms
        double epsilon; // kcal mol^-1
    
    public:
        // TODO Use eigen or armadillo for coordinates, but we will worry about that later
        double x;
        double y;
        double z;

        // af = analytical force
        double x_af = 0.0;
        double y_af = 0.0;
        double z_af = 0.0;

        Atom(int atomic_number_, double x_, double y_, double z_) :
            atomic_number(atomic_number_), x(x_), y(y_), z(z_)
            {
                if (atomic_number_ != 79) // Error handles non-gold. TODO handle < 1 cases
                {
                    throw std::exception();
                }
                sigma = 2.951;
                epsilon = 5.29;
            }
        
        int get_atomic_number() const { return atomic_number; }

        // TODO Change to armadillo vector
        void get_coords() const {}

        double get_sigma() const { return sigma; }
        double get_epsilon() const { return epsilon; }

        void set_atomic_number(int new_atomic_number)
        {
            atomic_number = new_atomic_number;
        }

        void set_coords()
        {

        }

        void const print_af() const
        {
            std::cout << x_af << " " << y_af << " " << z_af << std::endl;
        }

        // Atom operator-(Atom atom_2)
        // {
        //     x - atom_2.x;
        // }
};

// TODO Getter and Setter Functions
class Cluster
{
    private:
        std::vector<Atom> atoms;
        // Maybe use a vector for forces or maybe a matrix

    public:
        // TODO Change to allow file to be read through command line
        // TODO Allow for an entire folder to be run too
        bool load_atoms(std::string file)
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
        void add_atom(int atomic_number, double x, double y, double z)
        {
            Atom new_atom(atomic_number, x, y, z);
            atoms.push_back(new_atom);
        }

        int get_num_atoms() const { return atoms.size(); }

        std::vector<Atom> get_atoms() const { return atoms; }
        std::vector<Atom>& get_atoms() { return atoms; }

        void print_atoms()
        {
            for (auto& atom : atoms)
            {
                std::cout << atom.get_atomic_number() << "(" << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
            }
        }

        void print_analytical_force()
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
};

// Load atoms in a system
// Check boundaries?

// ENERGY

class Calculator
{
    public:
            // Utilizing mixing rules
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

};

class Energy_calc : public Calculator
{
    public:
        double calculate_energy(Cluster clusters)
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
};


// FORCE -> F = derivative of E / derivative of R?
// So take the derivative of LJ_potential for analytical force
// Use the function with step size difference (this means changes in radius?)

class Force_calc : public Calculator
{
    public: 
        // TODO: Worry about this later but this isn't working
        double calculate_analytical_force(Cluster& clusters)
        {
            int num_atoms = clusters.get_num_atoms();
            
            double total_force = 0;
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

            return total_force;
        }

        double forward_difference(double step_size)
        {
            // [f(x + h) - f(x)] / h

            return 1;
        }

        double central_difference(double step_size)
        {
            return 1;
        }

        double approximate_force(std::string finite_difference_method, double step_size)
        {
            if (finite_difference_method == "forward")
            {
                return forward_difference(step_size);
            }
            
            if (finite_difference_method == "central")
            {
                return central_difference(step_size);
            }
        }

};




int main(void)
{
    // 1. Create cluster
    Cluster gold;
    gold.load_atoms("./sample_input/Force/1.txt");

    // 2. Calculate Energy
    Energy_calc e_calc;
    gold.print_atoms();
    double energy_of_the_cluster = e_calc.calculate_energy(gold); // calculate_energy();

    // 3. Calculate Force
    Force_calc f_calc;
    double step_size = 0.1;
    double forw_diff = f_calc.forward_difference(step_size);
    double cent_diff = f_calc.central_difference(step_size);
    double analytical_force_of_the_cluster = f_calc.calculate_analytical_force(gold);
    gold.print_analytical_force();


    int force_of_the_cluster = 0;
    std::cout << "" << std::endl;

    std::cout << "Done" << std::endl;
    return 0;
}