#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "atom.hpp"
#include "lennard_jones.hpp"
#include <vector>
#include <string>

class Cluster
{
    private:
        std::vector<Atom> atoms;
        // Maybe use a vector for forces or maybe a matrix

    public:
        // TODO Change to allow file to be read through command line
        // TODO Allow for an entire folder to be run too
        bool load_atoms(std::string file);
        // Cluster(Cluster new_cluster, double step_size);
        // TODO: Find out if I have to template this or if I can just take in an integer, float, double etc.
        void add_atom(int atomic_number, double x, double y, double z);
        int get_num_atoms() const;
        std::vector<Atom> get_atoms() const;
        std::vector<Atom>& get_atoms();
        void print_atoms();
        double calculate_total_energy();
        void update_analytical_force();
        void print_analytical_force();
        void print_forward_difference();

};

#endif