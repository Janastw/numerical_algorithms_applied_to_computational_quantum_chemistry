#ifndef DERIVATIVE_APPROXIMATION_HPP
#define DERIVATIVE_APPROXIMATION_HPP

double forward_difference(double f_x, double f_x_plus_step_size, double step_size);

double central_difference(double f_x_plus_step_size, double f_x_minus_step_size, double step_size);


#endif