#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "atom.hpp"
#include "lennard_jones.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <functional>

class Cluster
{
    private:
        std::vector<Atom> atoms;
        // Maybe use a vector for forces or maybe a matrix

    public:
        arma::mat system_coordinates;
        arma::mat system_coordinate_forces;
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
        void print_system_coordinates();
        void print_system_coordinate_forces();
        double calculate_total_energy(const arma::mat& coords);
        void zero_forces();
        void update_analytical_force();
        void step_size(double step_size);
        void update_forward_difference(double step_size);
        void update_central_difference(double step_size);

        arma::vec gradient(int atom1_index, int atom2_index, double step_size);

        void steepest_descent(double step_size, double sd_step_size, double threshold);
        void print_analytical_force();
        double calculate_total_energy_at(double x);
        std::tuple<double, double, double> bracketing(double b, double step_size);
        void bracket(double a, double b, std::function<double(double)> operation, double& ax, double& bx, double& cx);
        void print_forward_difference();
        void print_steepest_descent();

};

struct Golden {
    double ax, bx, cx;
    double xmin, f_min;
    const double tol;

    Golden(double tol = 0.3) : tol(tol) {}

    template <class T>
    double minimize(T &func);

private:
    void shift3(double &a, double &b, double &c, double d) {
        a = b;
        b = c;
        c = d;
    }

    void shift2(double &a, double &b, double c) {
        a = b;
        b = c;
    }
};

#endif