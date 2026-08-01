#ifndef GMOM_H
#define GMOM_H

/* gmom - Generalized Method of Moments (divide-and-conquer), b=1 branch.
 *
 * Port of the reference MATLAB mbmom1.m / setup1.m (Casale, "A Generalized
 * Method of Moments for Closed Queueing Networks", Perf. Eval. 2011).
 *
 * The original MoM (bin/mom) recurses on the population only.  The
 * generalized MoM adds the Convolution recursion, which recurses on the
 * NUMBER OF QUEUES: a model on the prefix {1..m} of queues is built from
 * the model on {1..m-1} plus a single queue-removal (subtractive
 * convolution) branch.  This is the b=1 case: one queue is removed at a
 * time, along the prefix chain 1 -> 2 -> ... -> M.
 *
 * The per-level linear system is overdetermined; it is solved in exact
 * rational arithmetic by normal equations (A^T A) x = A^T b, i.e. eq (8)
 * of the paper.  gmpla supplies the sparse/dense mpq machinery and util
 * the combinatorial helpers, so this file reuses the same infrastructure
 * as bin/mom and bin/comom.
 */

#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include <gmpla.h>
#include "util.h"
#include "pfqn_sing.h"   /* LS1, setup1/ls1_free, singularity oracle (shared) */

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
extern struct rusage ruse;

/* LS1, setup1(), ls1_free() and the singularity oracle now live in
 * util/pfqn_sing.{c,h}: they are loading-only and shared with the rest of
 * the MoM family (see the scope note in that header). */

/* momblvl: the level indices used by the reference.  For b=1 the prefix
 * chain uses lvl=r-1 at slot l and lvl=r-2 at slot l_1 (clamped at 0). */
static inline int gmom_lvl(int r)   { return r-1; }
static inline int gmom_lvl_1(int r) { int v=r-2; return v<0?0:v; }

/* Capture mode for gmom_measures.  When a gmom_out* is passed, the measures
 * are written there as doubles instead of printed; -b/--bounds uses this to
 * bracket two perturbed solves that differ only in the perturbation seed.
 * X is indexed by original class, Q row-major as Q[(k-1)*R + r]. */
typedef struct {
	double  G;
	double  logG;
	double* X;   /* R entries    */
	double* Q;   /* M*R entries  */
} gmom_out;

/* base normalizing constant for the single-queue prefix (m=1): a model of
 * `mult` identical queues with demand row d[0..r-1] plus think times Z,
 * at population pop[0..r-1].  Computed exactly by population-lattice
 * convolution, valid for Z=0 and Z!=0 alike.  Result in G (mpq). */
void gmom_base(mpq_t G, mpz_t* d, int mult, int* pop, mpz_t* Z, int r);

#endif
