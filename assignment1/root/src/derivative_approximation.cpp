#include "derivative_approximation.hpp"

double forward_difference(double f_x, double f_x_plus_step_size, double step_size)
{
    // [f(x + h) - f(x)] / h
    double diff = f_x_plus_step_size - f_x;
    double division = diff / step_size;

    return division;
}

double central_difference(double f_x_plus_step_size, double f_x_minus_step_size, double step_size)
{
    // [f(x + h) - f(x - h)] / 2h
    return (f_x_plus_step_size - f_x_minus_step_size) / (2 * step_size);
}