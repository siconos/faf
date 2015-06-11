#ifndef BOGUS_INTERFACE_HPP
#define BOGUS_INTERFACE_HPP

#include <SolverOptions.h>
#include <Interfaces/FrictionProblem.hpp>

extern "C"
{
#include <fclib.h>
}

double norm_fb(double *reaction, double *velocity, double mu);

void compute_fb(double *reaction, double *velocity, double mu, double *out3);

void compute_jfb(double *reaction, double *velocity, double mu,
                 double *out3, double *out9_1, double *out9_2);

int solve_fclib( const fclib_local* problem, double* reactions, double* velocities, SolverOptions* SO );


#endif
