#ifndef FUNCODEGEN
#define FUNCODEGEN

/* Definitions are incomplete in frama-c Magnesium, so we use our
 * specs for math functions */
#ifndef __FRAMAC__
#include <math.h>
#endif

#include <float.h>
#include <assert.h>

/*@
  axiomatic Max {
  logic real fmax(real x, real y);
  axiom fmax_1: \forall real x, real y; fmax(x, y) >= x && fmax(y, y) >= y;
  axiom fmax_2: \forall real x, real y; fmax(x, y) == x || fmax(y, y) == y;
  }
  axiomatic sqrt {
  logic real sqrt(real x);
  axiom sqrt_1: \forall real x; x >= 0 <==> x == sqrt(x) * sqrt(x);
  axiom sqrt_2: \forall real x; x > 0 <==> sqrt(x) > 0;
  }
  axiomatic Heaviside {
  logic real Heaviside(real x);
  axiom Heaviside_1: \forall real x; (x < 0 ==> Heaviside(x) == 0);
  axiom Heaviside_2: \forall real x; (x > 0 ==> Heaviside(x) == 1);
  axiom Heaviside_3: \forall real x; (x == 0 ==> 0 < Heaviside(x) < 1);
  }
  axiomatic pow {
  logic real pow(real x, real y);
  axiom pow_1: \forall real x, real y; x >=0 ==> pow(x, y) >= 0;
  }
  axiomatic general {
  axiom sq: \forall real x; x*x >= 0.;
  }
*/

#define Sign(x) ((double)(x>0) - (double)(x<0))
#define Max fmax
#define Heaviside(x) (x < 0 ? 0 : ((x > 0) ? 1 : .5))
#define Rand(x) ((double) rand()/ (double) RAND_MAX)

/*@ requires \is_finite((double) x) && x >= 0.;
    assigns \nothing;
    ensures \is_finite((double) \result);
    ensures 0. <= \result <= 1. + x;
    ensures x <= \result <= 1. || 1. <= \result <= x;
    ensures \forall real a; \result >= a ==> x >= a*a;
    ensures \forall real a;  a <= x <= 1. ==> a*a;
    ensures \result * \result == x; */
extern double sqrt(double x);

/*@ requires \is_finite((double) x) && \is_finite((double) y);
    assigns \nothing;
    ensures \is_finite(\result) && (\result == x || \result == y) && (\result >= x && \result >= y);
 */
extern double fmax(double x, double y);

#ifdef __FRAMAC__
#include <__fc_builtin.h>
#endif

#endif
