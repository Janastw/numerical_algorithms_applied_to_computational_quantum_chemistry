#ifndef ATOM_HPP
#define ATOM_HPP

#include <iostream>
#include <string>
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
        double x, y, z;

        // af = analytical force
        double x_af = 0.0, y_af = 0.0, z_af = 0.0;

        Atom(int atomic_number_, double x_, double y_, double z_);

        int get_atomic_number() const;
        double get_sigma() const;
        double get_epsilon() const;

};

#endif