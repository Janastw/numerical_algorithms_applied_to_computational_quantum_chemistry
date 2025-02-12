#ifndef STEEPEST_DESCENT_HPP
#define STEEPEST_DESCENT_HPP

#include "derivative_approximation.hpp"
#include "cluster.hpp"
#include "lennard_jones.hpp"
#include <armadillo>

void steepest_descent_central_diff(Cluster cluster, double derivative_step_size, double sd_step_size, double convergence_threshold, int iteration);


#endif