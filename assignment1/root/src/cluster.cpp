#include "cluster.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <armadillo>



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

    arma::mat system_coordinates(3, num_atoms, arma::fill::zeros);

    int count = 0;

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
    system_coordinate_forces = arma::mat(3, num_atoms, arma::fill::zeros);
    return true;
}

void Cluster::add_atom(int atomic_number, double x, double y, double z)
{
    Atom new_atom(atomic_number, x, y, z);
    atoms.push_back(new_atom);
    system_coordinates = arma::join_rows(system_coordinates, new_atom.coords);
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

void Cluster::print_system_coordinates()
{
    std::cout << std::fixed << std::setprecision(10) << system_coordinates << std::endl;
}

void Cluster::print_system_coordinate_forces()
{
    std::cout << std::fixed << std::setprecision(10) << system_coordinate_forces << std::endl;
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
    for (int i = 0; i < 3; i++)
    {
        for (auto& atom: atoms)
        {
            std::cout << std::setw(8) << std::setprecision(4) << atom.coords_analytical_forces[i] << " ";
        }
        std::cout << std::endl;
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
            double sigma_ij = calculate_sigma_ij(atoms[i].get_sigma(), atoms[j].get_sigma());
            double radius_ij = calculate_distance(system_coordinates.col(i), system_coordinates.col(j));
            double epsilon_ij = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[j].get_epsilon());
            total_energy += calculate_lj_energy(sigma_ij, radius_ij, epsilon_ij);
        }
    }

    return total_energy;
}

void Cluster::zero_forces()
{
    for (auto& atom: atoms)
    {
        atom.coords_analytical_forces = arma::vec({0.0, 0.0, 0.0});
    }
    system_coordinate_forces.zeros();
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
            double radius_ik = calculate_distance(atoms[i].coords, atoms[k].coords);
            double epsilon_ik = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[k].get_epsilon());
            double analytical_force = calculate_lj_force(sigma_ik, radius_ik, epsilon_ik);

            arma::vec direction = system_coordinates.col(i) - system_coordinates.col(k);
            system_coordinate_forces.col(i) += analytical_force * direction;

            // Can remove this line when my matrix forces are well integrated
            atoms[i].coords_analytical_forces += (atoms[i].coords - atoms[k].coords) * analytical_force;
        }
    }
}

void Cluster::update_forward_difference(double step_size)
{
    int num_atoms = atoms.size();
    std::vector<Atom>& atoms = get_atoms();

    for (int i = 0; i < num_atoms; i++)
    {
        for (int j = 0; j < num_atoms; j++)
        {
            if (i == j)
            {
                continue;
            }
            double sigma_ij = calculate_sigma_ij(atoms[i].get_sigma(), atoms[j].get_sigma());
            double epsilon_ij = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[j].get_epsilon());
            double radius_ij = calculate_distance(atoms[i].coords, atoms[j].coords);

            double f_of_x = calculate_lj_energy(sigma_ij, radius_ij, epsilon_ij);
            double f_of_x_plus_step_size = calculate_lj_energy(sigma_ij, radius_ij + step_size, epsilon_ij);

            double forward_diff = forward_difference(f_of_x, f_of_x_plus_step_size, step_size);
            
            arma::vec direction = arma::normalise(system_coordinates.col(i) - system_coordinates.col(j));
            system_coordinate_forces.col(i) += forward_diff * direction;

            // Can remove this line when my matrix forces are well integrated
            atoms[i].coords_analytical_forces += arma::normalise(atoms[i].coords - atoms[j].coords) * forward_diff;

        }
    }
}

void Cluster::update_central_difference(double step_size)
{
    int num_atoms = atoms.size();
    std::vector<Atom>& atoms = get_atoms();

    for (int i = 0; i < num_atoms; i++)
    {
        for (int j = 0; j < num_atoms; j++)
        {
            if (i == j)
            {
                continue;
            }
            double sigma_ij = calculate_sigma_ij(atoms[i].get_sigma(), atoms[j].get_sigma());
            double epsilon_ij = calculate_epsilon_ij(atoms[i].get_epsilon(), atoms[j].get_epsilon());
            double radius_ij = calculate_distance(atoms[i].coords, atoms[j].coords);

            double f_of_x_plus_step_size = calculate_lj_energy(sigma_ij, radius_ij + step_size, epsilon_ij);
            double f_of_x_minus_step_size = calculate_lj_energy(sigma_ij, radius_ij - step_size, epsilon_ij);
            double central_diff = central_difference(f_of_x_minus_step_size, f_of_x_plus_step_size, step_size);

            arma::vec direction = arma::normalise(system_coordinates.col(i) - system_coordinates.col(j));
            system_coordinate_forces.col(i) += central_diff * direction;

            // Can remove this line when my matrix forces are well integrated
            atoms[i].coords_analytical_forces += arma::normalise(atoms[i].coords - atoms[j].coords) * central_diff;
        }
    }
}

template <class T>
double Golden::minimize(T &func) {
    const double R = 0.61803399;  // Golden ratio
    const double C = 1.0 - R; 

    double x0 = ax, x3 = cx;
    double x1, x2;

    // Initialize x1 and x2 within the bracket
    if (fabs(cx - bx) > fabs(bx - ax)) {
        x1 = bx;
        x2 = bx + C * (cx - bx);
    } else {
        x2 = bx;
        x1 = bx - C * (bx - ax);
    }
    double f1 = func(x1);
    double f2 = func(x2);

    while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
        if (f2 < f1) {
            shift3(x0, x1, x2, R * x2 + C * x3);
            shift2(f1, f2, func(x2));
        } else {
            shift3(x3, x2, x1, R * x1 + C * x0);
            shift2(f2, f1, func(x1));
        }
    }

    if (f1 < f2) {
        xmin = x1;
        f_min = f1;
    } else {
        xmin = x2;
        f_min = f2;
    }

    return xmin;
}

arma::vec Cluster::gradient(int atom1_index, int atom2_index, double step_size)
{
    double sigma_ij = calculate_sigma_ij(atoms[atom1_index].get_sigma(), atoms[atom2_index].get_sigma());
    double epsilon_ij = calculate_epsilon_ij(atoms[atom1_index].get_epsilon(), atoms[atom2_index].get_epsilon());
    double radius_ij = calculate_distance(system_coordinates.col(atom1_index), system_coordinates.col(atom2_index));

    double f_of_x_plus_step_size = calculate_lj_energy(sigma_ij, radius_ij + step_size, epsilon_ij);
    double f_of_x_minus_step_size = calculate_lj_energy(sigma_ij, radius_ij - step_size, epsilon_ij);

    double scalar = -1 * central_difference(f_of_x_plus_step_size, f_of_x_minus_step_size, step_size);

    arma::vec direction = arma::normalise(system_coordinates.col(atom1_index) - system_coordinates.col(atom2_index));
    // system_coordinate_forces.col(i) += central_diff * direction;

    arma::vec force = scalar * direction;

    return force;
}

void Cluster::steepest_descent(double step_size, double sd_step_size, double threshold)
{
    arma::mat new_coords = system_coordinates;
    double old_energy = calculate_total_energy();
    for (int i = 0; i < system_coordinates.n_cols; i++)
    {
        arma::vec total_gradient = arma::zeros<arma::vec>(system_coordinates.n_rows);
        for (int j = 0; j < system_coordinates.n_cols; j++)
        {
            if (i == j)
            {
                continue;
            }
            total_gradient += gradient(i, j, step_size);
        }

        // auto [a, b, c] = bracketing(sd_step_size, threshold);

        // Golden golden; golden.ax = a; golden.bx = b; golden.cx = c;
        // auto energy_function = [this](double x) { return calculate_total_energy_at(x);};
        // double new_step_size = golden.minimize(energy_function);
        double a, b, c;
        bracket(sd_step_size, threshold, [this](double x) { return calculate_total_energy_at(x); }, a, b, c);
        Golden golden; 
        golden.ax = a; 
        golden.bx = b; 
        golden.cx = c;
        auto energy_function = [this](double x) { return calculate_total_energy_at(x); };
        double new_step_size = golden.minimize(energy_function);


        new_coords.col(i) -= new_step_size * total_gradient;
    }

    system_coordinates = new_coords;
    update_central_difference(step_size);

    double new_energy = calculate_total_energy();

    if (std::abs(new_energy - old_energy) < threshold)
    {
        return;
    }
}

double Cluster::calculate_total_energy_at(double x)
{
    arma::mat temp_coordinates = system_coordinates;

    for (int i = 0; i < temp_coordinates.n_cols; i++)
    {
        arma::vec total_gradient = arma::zeros<arma::vec>(temp_coordinates.n_rows);
        for (int j = 0; j < temp_coordinates.n_cols; j++)
        {
            if (i == j) continue;
            total_gradient += gradient(i, j, x);
        }

        temp_coordinates.col(i) -= x * total_gradient;
    }

    return calculate_total_energy();
}

std::tuple<double, double, double> Cluster::bracketing(double b, double step_size)
{
    double a = b - step_size;
    double c = b + step_size;

    double f_a = calculate_total_energy_at(a);
    double f_b = calculate_total_energy_at(b);
    double f_c = calculate_total_energy_at(c);

    while (!(f_b < f_a && f_b < f_c))
    {
        if (f_a < f_b)
        {
            c = b;
            b = a;
            a -= step_size;
            f_c = f_b;
            f_b = f_a;
            f_a = calculate_total_energy_at(a);
        }
        else if (f_c < f_b)
        {
            a = b;
            b = c;
            c += step_size;
            f_a = f_b;
            f_b = f_c;
            f_c = calculate_total_energy_at(c);
        }
        else
        {
            break;
        }
    }

    return {a, b, c};
}

void Cluster::bracket(double a, double b, std::function<double(double)> operation, double& ax, double& bx, double& cx)
{
    const double GOLD = 1.618034;  // Golden ratio constant
    const double GLIMIT = 100.0;
    const double TINY = 1.0E-20;

    ax = a;
    bx = b;
    double fa = operation(ax);
    double fb = operation(bx);

    // Ensure ordering such that f(a) > f(b)
    if (fb > fa)
    {
        std::swap(ax, bx);
        std::swap(fa, fb);
    }

    // First guess for c
    cx = bx + GOLD * (bx - ax);
    double fc = operation(cx);

    while (fb > fc)  // Keep expanding the bracket
    {
        double r = (bx - ax) * (fb - fc);
        double q = (bx - cx) * (fb - fa);
        double denom = 2.0 * std::max(std::abs(q - r), TINY);
        double u = bx - ((bx - cx) * q - (bx - ax) * r) / denom;
        double ulim = bx + GLIMIT * (cx - bx);
        double fu;

        if ((bx - u) * (u - cx) > 0.0)  // Parabolic fit is between b and c
        {
            fu = operation(u);
            if (fu < fc)
            {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            }
            else if (fu > fb)
            {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + GOLD * (cx - bx);
            fu = operation(u);
        }
        else if ((u - ulim) * (ulim - cx) >= 0.0)  // Limit parabolic fit to max bracket size
        {
            u = ulim;
            fu = operation(u);
        }
        else  // Default golden ratio step
        {
            u = cx + GOLD * (cx - bx);
            fu = operation(u);
        }

        // Shift the bracket
        ax = bx;
        bx = cx;
        cx = u;
        fa = fb;
        fb = fc;
        fc = fu;
    }
}

