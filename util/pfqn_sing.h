#ifndef PFQN_SING_H
#define PFQN_SING_H

/* Shared level construction and singularity oracle for the MoM family.
 *
 * setup1() builds the b=1 generalized-MoM level system; its coefficient
 * matrix depends on the loadings L only, not on N or Z.  pfqn_sing.c uses
 * that property to decide, once at initialisation, whether the recursion
 * would hit a rank-deficient matrix, and which class to move to the
 * recursion position to avoid it.
 *
 * Scope note: the oracle is exact for the b=1 matrix built here (gmom).
 * bin/mom and bin/comom use a different basis, so for those it is a hint
 * about the underlying loading degeneracy, not a decision procedure --
 * candidate class orders must still be probed by running the solver
 * (see util/safe_tiers.c).
 */

#include <gmp.h>
#include <gmpla.h>
#include "util.h"

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

/* singularity oracle (pfqn_sing.c) */
int pfqn_recursion_singular(mpz_t** L, int M, int R);
int pfqn_nonsingular_recclass(mpz_t** L, int M, int R);

#endif
