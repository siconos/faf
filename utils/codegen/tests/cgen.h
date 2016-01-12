//#include <math.h>
#include <float.h>
#include <assert.h>

/*@
  axiomatic Max {
  logic real fmax(real x, real y);
  axiom fmax_1: \forall real x, real y; fmax(x, y) >= x && fmax(y, y) >= y;
  axiom fmax_2: \forall real x, real y; fmax(x, y) == x || fmax(y, y) == y;
  }
  axiomatic Sqrt {
  logic real sqrt(real x);
  axiom sqrt_1: \forall real x; x >= 0 <==> x == sqrt(x) * sqrt(x);
  axiom sqrt_2: \forall real x; x > 0 <==> sqrt(x) > 0;
  }
  axiomatic Mul {
  logic real mul(real x, real y);
  axiom mul_1: \forall real x, real y; x>=0 && y>=0 ==> x * y >= 0;
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
  axiom div: \forall real x; x != 0 ==> \is_finite((double) (1/x));
  }

*/

/*@ predicate same_sign(real x, real y) = ((x > 0) <==> (y > 0)) && ((x < 0) <==> (y < 0)) && ((x == 0) <==> (y == 0));
 */
#define DEBUG_PRINT(X)
#define Sign(x) ((double)(x>0) - (double)(x<0))
#define Max fmax
#define Heaviside(x) (x < 0 ? 0 : ((x > 0) ? 1 : .5))
#define Rand(x) ((double) rand()/ (double) RAND_MAX)

/* Round-off errors: for values supposed to be >=0,  value = 0 is replaced by value <= ZERO */
/* Note:
 * ZERO=DBL_EPSILON*100 (too small) => test-3D_NSN_FB-NESpheres_30_1 fails (division by 0)
 * ZERO=1e-10 (too big) => test-3D_NSN_FB-OneObject-i100000-499.hdf5 fails (no convergence)
 */
#define ZERO DBL_EPSILON*500

#define POST_CHECK_POW(x)
#define POST_CHECK_ADD(x)
#define POST_CHECK_MUL(x)
#define POST_CHECK(x)

#define NOT_ZERO(x) fabs(x) > 0
#define IS_NOT_ZERO(x) fabs(x) > 0
#define IS_POSITIVE(x) 1



#define random1 sqrt(2)/2
#define random2 sqrt(2)/2

#ifdef __cplusplus
#include <cmath>
#define CHECK(x)
#define XCHECK(x) assert(isfinite(x))
#else
#define CHECK(x)
#define XCHECK(x) assert(isfinite(x))
#endif

#pragma GCC diagnostic ignored "-Wconversion"

// hack, should be prevented in sage/sympy/maple or in code generation
//#define sqrt(x) ((fabs((double)x) <= ZERO) ? 0 : (assert(x>=0), sqrt(x)))


/*@ requires \is_finite(x) && x >= 0;
    assigns \nothing;
    ensures \is_finite(\result)  && (0 <= \result <= x) && (x > 0 <==> \result > 0);
 */
extern double sqrt(double x);

/*@ requires \is_finite(x) && \is_finite(y);
    assigns \nothing;
    ensures \is_finite(\result);
 */
extern double mul(double x, double y);


/*@ requires \is_finite(x) && \is_finite(y);
    assigns \nothing;
    ensures \is_finite(\result) && (\result == x || \result == y) && (\result >= x && \result >= y);
 */
extern double fmax(double x, double y);





/*@ requires \is_finite(x) && \is_finite(y);
    assigns \nothing;
    ensures \is_finite(\result) && (x >= 0 ==> \result >= 0);
 */
extern double pow(double x, double y);

#include <builtin.h>
double Frama_C_double_interval(double, double);

