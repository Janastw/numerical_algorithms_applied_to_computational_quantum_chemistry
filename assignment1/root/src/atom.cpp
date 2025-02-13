#include "atom.hpp"

Atom::Atom(int atomic_number_, double x_, double y_, double z_) :
    atomic_number(atomic_number_), x(x_), y(y_), z(z_)
    {
        if (atomic_number_ != 79)
        {
            throw std::exception();
        }
        sigma = 2.951;
        epsilon = 5.29;
        coords = {x_, y_, z_};
    }

int Atom::get_atomic_number() const { return atomic_number; }
double Atom::get_sigma() const { return sigma; }
double Atom::get_epsilon() const { return epsilon; }