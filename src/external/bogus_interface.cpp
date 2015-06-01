
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "bogus_interface.hpp"

#include <Core/Block.impl.hpp>
#include <Core/Block.io.hpp>
#include <Core/BlockSolvers/ProjectedGradient.impl.hpp>
#include <Core/BlockSolvers/GaussSeidel.impl.hpp>
#include <Interfaces/FrictionProblem.hpp>
#include <Interfaces/FrictionProblem.impl.hpp>
#include <Extra/SecondOrder.fwd.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>


/* This is for a switch over bogus strategies : PureNewton,
 * PureEnumerative, Hybrid, RevHybrid (default for coulomb
 * friction). */
namespace bogus {
  template< unsigned Dimension, bogus::local_soc_solver::Strategy Strat >
  struct DualFrictionProblemInterf;

  namespace friction_problem {
    template< unsigned Dimension, template <typename> class Method, bogus::local_soc_solver::Strategy Strat >
    static double xsolve( const DualFrictionProblemInterf< Dimension, Strat >& dual,
                          const ConstrainedSolverBase< Method, typename DualFrictionProblemInterf< Dimension, Strat>::WType > &gs,
                          double *r, const bool staticProblem );
  }
  template< unsigned Dimension, bogus::local_soc_solver::Strategy Strat >
  struct DualFrictionProblemInterf : public DualFrictionProblem<Dimension>
  {
    typedef SOCLaw< Dimension, double, true, Strat  > CoulombLawType	;
    typedef SOCLaw< Dimension, double, false, Strat > SOCLawType	;

    double xsolveWith( typename DualFrictionProblem<Dimension>::GaussSeidelType &gs, double * r, const bool staticProblem = false ) const
    {
      gs.setMatrix( this->W );
      return bogus::friction_problem::xsolve( *this, gs, r, staticProblem ) ;
    }
  };


  namespace friction_problem {

    template< unsigned Dimension, template <typename> class Method, bogus::local_soc_solver::Strategy Strat >
    static double xsolve( const DualFrictionProblemInterf< Dimension, Strat >& dual,
                          const ConstrainedSolverBase< Method, typename DualFrictionProblemInterf< Dimension, Strat>::WType > &gs,
                          double *r, const bool staticProblem )
    {
      typename Eigen::VectorXd::MapType r_map ( r, dual.W.rows() ) ;

      if( dual.permuted() )
        applyPermutation< Dimension >( dual.permutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

      double res = staticProblem
        ? gs.solve( typename DualFrictionProblemInterf< Dimension, Strat >::SOCLawType
                    ( dual.W.rowsOfBlocks(), dual.mu.data() ), dual.b, r_map )
        : gs.solve( typename DualFrictionProblemInterf< Dimension, Strat >::CoulombLawType
                    ( dual.W.rowsOfBlocks(), dual.mu.data() ), dual.b, r_map ) ;

      if( dual.permuted() )
        applyPermutation< Dimension >( dual.invPermutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

      return res ;
    }

  }


}
//#define USE_PG_FOR_CADOUX


static const char* g_meth  = 0;

unsigned int g_gs_iter;
double g_gs_err;

void ackCurrentResidual( unsigned GSIter, double err )
{
  g_gs_iter = GSIter;
  g_gs_err = err;
//	std::cout << " .. " << g_meth << ": " <<  GSIter << " ==> " << err << std::endl ;
}

template< unsigned Dimension, bogus::local_soc_solver::Strategy Strat, typename EigenDerived >
static double solve( const fclib_local* problem, const Eigen::SparseMatrixBase< EigenDerived >& ei_W,
                     Eigen::VectorXd &r, Eigen::VectorXd &u, const bool useCadoux, SolverOptions* SO )
{
  double tol = SO->dparam[0];

	bogus::DualFrictionProblemInterf< Dimension, Strat > dual ;
	bogus::convert( ei_W, dual.W, Dimension, Dimension ) ;

	dual.W.prune( tol ) ;
	dual.W.cacheTranspose();

	dual.b = Eigen::VectorXd::Map( problem->q, problem->W->n ) ;
	dual.mu = Eigen::VectorXd::Map( problem->mu, problem->W->n/Dimension ) ;

	typename bogus::DualFrictionProblem< Dimension >::GaussSeidelType gs ;
	gs.setTol( tol ) ;
  gs.setEvalEvery( 1 ) ;
	gs.setAutoRegularization( 1.e-3 ) ;

	bogus::Signal< unsigned, double > callback ;
	callback.connect( &ackCurrentResidual );

	double res = -1 ;

	if( useCadoux )
	{
		g_meth = "Cadoux" ;
#ifdef USE_PG_FOR_CADOUX
		typename bogus::DualFrictionProblem< Dimension >::ProjectedGradientType pg ;
		pg.setTol( tol ) ;
		pg.setMaxIters( 1000 ) ;
		res = dual.solveCadoux( pg, r.data(), 100, &callback ) ;
#else
		gs.setMaxIters( 1000 ) ;
		res = dual.solveCadoux( gs, r.data(), 100, &callback ) ;
#endif

	} else {
		g_meth = "GS" ;
		gs.setMaxIters( 1.e5 );
		gs.callback().connect( callback );

 		res = dual.xsolveWith( gs, r.data() ) ;
	}

	u = dual.W * r + dual.b ;

  SO->dparam[1] = g_gs_err;
  SO->iparam[1] = g_gs_iter;
	return res ;
}

int solve_fclib( const fclib_local* problem, SolverOptions* SO )
{

  double tolerance = SO->dparam[0];

  int strategy = SO->iparam[4];

  bool useCadoux = false;

  if( problem )
  {
		/*std::cout << "Successfully loaded problem " << file << std::endl ;
		std::cout << " Name: " << problem->info->title << std::endl ;
		std::cout << " Descr: " << problem->info->description << std::endl ;
		std::cout << " Info: " << problem->info->math_info << std::endl ;

		if( problem->V ) std::cout << " Has V " << std::endl ;
		if( problem->R ) std::cout << " Has R " << std::endl ;*/

		const unsigned d = problem->spacedim  ;
		const unsigned n = problem->W->m / d  ;

		if( problem->W && !problem->V && !problem->R )
		{
			/*std::cout << " Pure " << d << "D Coulomb friction problem with "
        << n << " contacts " << std::endl ;*/

			if( problem->W->nz == -2 )
			{
				/*std::cout << " Compressed row storage " << problem->W->m << " / " << problem->W->n << " / " << problem->W->nzmax << std::endl ;*/

				Eigen::SparseMatrix< double, Eigen::RowMajor > ei_W ;
				ei_W.resize( problem->W->m, problem->W->n );
				ei_W.resizeNonZeros( problem->W->nzmax ) ;

				memcpy( ei_W.outerIndexPtr(), problem->W->p, ( problem->W->m+1 ) * sizeof( int ) ) ;
				memcpy( ei_W.innerIndexPtr(), problem->W->i, problem->W->nzmax * sizeof( int ) ) ;
				memcpy( ei_W.valuePtr(), problem->W->x, problem->W->nzmax * sizeof( double ) ) ;

				if( 0 == problem->W->p[ problem->W->m ] )
				{
					std::cout << " /!\\ Malformed sparse matrix ; entering repair mode " << std::endl ;
					assert( problem->W->nzmax == problem->W->n * problem->W->m ) ;

					ei_W.outerIndexPtr()[ 0 ] = 0 ;
					for( int row = 0 ; row < problem->W->m ; ++row )
					{
						const int start = ei_W.outerIndexPtr()[ row ] ;
						ei_W.outerIndexPtr()[ row+1 ] = start + problem->W->n ;
						for( int col = 0 ; col < problem->W->n ; ++col )
						{
							ei_W.innerIndexPtr()[ start + col ] = col ;
						}
					}

				}


				double res = -1. ;
				Eigen::VectorXd r, u ;
				r.setZero( problem->W->n ) ;

        switch (strategy)
        {
        case 0:

          if( problem->spacedim == 3 )
          {
            res = solve< 3u, bogus::local_soc_solver::PureNewton>( problem, ei_W, r, u, useCadoux, SO ) ;
          } else {
            res = solve< 2u, bogus::local_soc_solver::PureNewton>( problem, ei_W, r, u, useCadoux, SO ) ;
          }
          break;

        case 1:
          if( problem->spacedim == 3 )
          {
            res = solve< 3u, bogus::local_soc_solver::PureEnumerative>( problem, ei_W, r, u, useCadoux, SO ) ;
          } else {
            res = solve< 2u, bogus::local_soc_solver::PureEnumerative>( problem, ei_W, r, u, useCadoux, SO ) ;
          }
          break;

        case 2:
          if( problem->spacedim == 3 )
          {
            res = solve< 3u, bogus::local_soc_solver::Hybrid>( problem, ei_W, r, u, useCadoux, SO ) ;
          } else {
            res = solve< 2u, bogus::local_soc_solver::Hybrid>( problem, ei_W, r, u, useCadoux, SO ) ;
          }
          break;

        case 3:
          if( problem->spacedim == 3 )
          {
            res = solve< 3u, bogus::local_soc_solver::RevHybrid>( problem, ei_W, r, u, useCadoux, SO ) ;
          } else {
            res = solve< 2u, bogus::local_soc_solver::RevHybrid>( problem, ei_W, r, u, useCadoux, SO ) ;
          }
          break;

        }

        return res >= tolerance;

      }
    }
  }
}
