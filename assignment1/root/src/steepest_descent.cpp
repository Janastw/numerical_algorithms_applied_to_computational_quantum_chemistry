#include "steepest_descent.hpp"
#include <armadillo>
#include <tuple>

void steepest_descent_central_diff(Cluster cluster, double derivative_step_size, double sd_step_size, double convergence_threshold, int iteration)
{
    if (iteration == 0)
    {
        std::cout << "start steepest descent with golden section line search" << "\n" <<
        "Initial energy: " << cluster.calculate_total_energy() << "\n" <<
        "Stepsize for central difference is:" << derivative_step_size <<
        ";Initial stepsize for line search is:" << sd_step_size <<
        ";Threshold for convergence in force is:" << "" << "\n" <<
        "Central Difference Force" << "\n" << std::endl;
        cluster.update_central_difference(derivative_step_size);
        cluster.print_forward_difference();
    }
    std::cout << "Start steepest descent with golden section line search using central difference force" << "\n" <<
    "Start golden section search\nnew_point" << std::endl;
    cluster.print_forward_difference();
    double original_energy = cluster.calculate_total_energy();
    std::cout << "current energy" << original_energy << std::endl;


    cluster.update_central_difference(derivative_step_size);
    double new_energy = cluster.calculate_total_energy();
    std::cout << "Central Difference Force" << std::endl;
    cluster.print_forward_difference();

    if (new_energy == original_energy)
    {
        std::cout << "Total iterations: " << iteration << "\n" <<
        "Final energy: " << new_energy << "\n" <<
        "Optimized structure:" << std::endl;
        cluster.print_atoms();

        return;
    }
    
    if (new_energy > original_energy)
    {
        steepest_descent_central_diff(cluster, derivative_step_size, sd_step_size / 2, convergence_threshold, iteration + 1);
    }
    else
    {
        steepest_descent_central_diff(cluster, derivative_step_size, sd_step_size * 1.0, convergence_threshold, iteration + 1);
    }
}

// template<typename T>
// std::tuple<double, double, double> bracket(const double x_0, T& f, double step_size, double expansion_factor=2)
// {
//     double a = x_0;
//     double b = a + step_size;
//     double c = b + step_size;

//     double f_a = f(a);
//     double f_b = f(b);
//     double f_c = f(c);

//     while (!(f_b < f_a && f_b < f_c))
//     {
//         a = b;
//         b = c;
//         c = c + step_size;
//         step_size += expansion_factor;

//         f_a = f(a);
//         f_b = f(b);
//         f_c = f(c);
//     }

//     return {a, b, c};
 
// }

// Golden golden;
// golden.bracker(a, b, func);
// xmin = golden.minimize(func);

// tol = ...
// Golden golden(tol);

// golden.ax = ...; golden.bx = ...; golden.cx = ...;

// struct Golden : Bracketmethod
// {
//     double xmin, fmin;
//     const double tol;
//     Golden(const double toll=3.0e-8) : tol(toll) {}
//     template <template T>
//     double minimize(T &func)
//     {
//         const double R = 0.61803399, C = 1.0-R;
//         double x1, x2;
//         double x0 = ax;
//         double x3 = cx;
//         if (abs(cx-bx) > abs(bx-ax))
//         {
//             x1=bx;
//             x2=bx + C * (cx - bx);
//         } else {
//             x2 = bx;
//             x1 = bx - C * (bx - ax);
//         }
//         double f1 = func(x1);
//         double f2 = func(x2);
//         while (abs(x3 - x0) > tol * (abs(x1) + abs(x2)))
//         {
//             if (f2 < f1) 
//             {
//                 shft3(x0, x1, x2, R * x2 + C * x3);
//                 shft2(f1, f2, func(x2));
//             } else {
//                 shft3(x3, x2, x1, R*x1 + C * x0);
//                 shft2(f2, f1, func(x1));
//             }
//         }
//         if (f1 < f2)
//         {
//             xmin = x1;
//             fmin = f1;
//         } else {
//             xmin = x2;
//             fmin = f2;
//         }
//         return xmin;
//     }
// }


void steepest_descent_2()
{

}

// struct Bracketmethod
// {
//     double ax, bx, cx, fa, fb, fc;
//     void bracket(double a, double b, double (*operation)(double, double, double))
//     {
//         double gold = 1.618034, glimit = 100.0, tiny = 1.0E-20;
//         ax = a;
//         bx = b;
//         fa = operation(a);
//         fb = operation(b);

//         if (fb > fa)
//         {
//             double temp1 = ax;
//             ax = bx;
//             bx = temp1;

//             double temp2 = fa;
//             fa = fb;
//             fb = temp2;
//         }

//         cx = bx + gold * (bx - ax);
//         fc = operation(cx);

//         while (fb > fc)
//         {
//             double r = (bx - ax) * (fb - fc);
//             double q = (bx - cx) * (fb - fa);
//             double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * sign(Max(abs(q-r), tiny), q-r));
//             double ulim = bx + glimit * (cx - bx);
//             if ((bx - u) * (u - cx) > 0.0) 
//             {
//                 fu = operation(u);
//                 {
//                     ax = bx;
//                     bx = u;
//                     fa = fb;
//                     fb = fu;
//                     return;
//                 } else if (fu > fb) {
//                     cx = u;
//                     fc = fu;
//                     return;
//                 }
//                 u = cx + gold * (cx -bx);
//                 fu = operation(u);
//             }
//             u = cx * gold * (cx - bx);
//             fu = operation(u);
//             if (fu < fc) 
//             {
//                 shft3(bx, cx, u, u + gold * (u - cx));
//                 shft3(fb, fc, fu, operation(u));
//             }
//         } else if ((u -ulim) * (ulim - cx) >= 0.0) {
//             u = ulim;
//             fu = operation(u);
//         } else {
//             u = cx + gold * (cx - bx);
//             fu = operation(u);
//         }
//         shft3(ax, bx, cx, u);
//         shft3(fa, fb, fc, fu);
//     }
//     inline void shft2(double &a, double &b, const double c)
//     {
//         a = b;
//         b = c;
//     }
//     inline void shft3(double &a, double &b, double &c, const double d)
//     {
//         a = b;
//         b = c;
//         c = d;
//     }
//     inline void mov3(double &a, double &b, double &c, const double d, const double e, const double f)
//     {
//         a =d; b = e; c = f;
//     }
// };




void golden_section_search();

arma::vec gradient_();

// arma::vec gradient(double f_of_x_plus_h, double f_of_x_minus_h, arma::vec coords)
// {
//     return -central_difference(f_of_x_plus_h, f_of_x_minus_h) * coords;
// }

void steepest_descent()
{
    // arma::vec gradient = central_difference()
}