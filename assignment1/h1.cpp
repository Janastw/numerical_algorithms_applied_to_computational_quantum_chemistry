// #include <armadillo>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>



class Atom
{
    private:
        int atomic_number;
        // TODO Use eigen or armadillo for coordinates, but we will worry about that later
        int x;
        int y;
        int z;
    
    public:
        Atom(int atomic_number_, int x_, int y_, int z_) :
            atomic_number(atomic_number_), x(x_), y(y_), z(z_) {}
        
        int get_atomic_number() const { atomic_number; }

        // TODO Change to armadillo vector
        void get_coords() const {}

        void set_atomic_number()
        {

        }

        void set_coords()
        {

        }
};

// TODO Getter and Setter Functions
class Simulation
{
    private:
        std::vector<Atom> atoms;
        int num_atoms = 0;

    public:
        void add_atom(int atomic_number, int x, int y, int z)
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
                std::cout << atom.get_atomic_number() << std::endl;
            }
        }




};

// Load atoms in a system
// Check boundaries?


// TODO Change to allow file to be read through command line
    // TODO Allow for an entire folder to be run too
Simulation load_atoms(std::string file)
{
    std::ifstream inputFile(file);
    std::string line;
    getline(inputFile, line);

    // First line is the number of atoms
    int num_atoms = std::stoi(line);
    std::cout << num_atoms << std::endl;

    class Simulation sim;

    // Input files have atomic number, x, y, then z coordinate in each line
    while (getline(inputFile, line))
    {
        int atomic_number;
        int x;
        int y;
        int z;
        
        std::istringstream numbers(line);
        numbers >> atomic_number;
        numbers >> x;
        numbers >> y;
        numbers >> z;

        sim.add_atom(atomic_number, x, y, z);
        
    }

    return sim;
}


int main(void)
{
    Simulation gold = load_atoms("./sample_input/Energy/1.txt");
    gold.print_atoms();
    std::cout << "Done" << std::endl;
    return 0;
}