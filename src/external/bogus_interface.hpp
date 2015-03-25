#ifndef BOGUS_INTERFACE_HPP
#define BOGUS_INTERFACE_HPP

#include <SolverOptions.h>
#include <Interfaces/FrictionProblem.hpp>

extern "C"
{
#include <fclib.h>
}

int solve_fclib( const fclib_local* problem, SolverOptions* SO );


#endif
