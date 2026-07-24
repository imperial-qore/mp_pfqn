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

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
extern struct rusage ruse;

/* The six per-level matrices produced by setup1, all dense mpq.
 * A,B are (nrows x szIk*R); C is (nrows x szIi*R); D is (nrows x szI*R);
 * E is (szI*R x szIk*R); F is (szI*R x szI*R).  nrows is the overdetermined
 * row count of the CE+PC system and is computed, not assumed. */
typedef struct
{
	mpq_mat_t A, B, C, D, E, F;
	int nrows;     /* rows of A,B,C,D (overdetermined) */
	int szIk, szIi, szI; /* combination counts at this level */
	int colA;      /* = szIk*R  (also rows of E,F) */
	int colC;      /* = szIi*R */
	int colD;      /* = szI*R */
} LS1;

/* setup1: build the six matrices for the prefix sub-model on queues 1..m
 * (queue m is the one removed to reach level m-1).  L is m x r, N and Z
 * are length r.  Returns a freshly allocated LS1. */
LS1* setup1(mpz_t** L, int m, int r, int* N, mpz_t* Z);
void ls1_free(LS1* ls, int r);

/* momblvl: the level indices used by the reference.  For b=1 the prefix
 * chain uses lvl=r-1 at slot l and lvl=r-2 at slot l_1 (clamped at 0). */
static inline int gmom_lvl(int r)   { return r-1; }
static inline int gmom_lvl_1(int r) { int v=r-2; return v<0?0:v; }

/* base normalizing constant for the single-queue prefix (m=1): a model of
 * `mult` identical queues with demand row d[0..r-1] plus think times Z,
 * at population pop[0..r-1].  Computed exactly by population-lattice
 * convolution, valid for Z=0 and Z!=0 alike.  Result in G (mpq). */
void gmom_base(mpq_t G, mpz_t* d, int mult, int* pop, mpz_t* Z, int r);

#endif
