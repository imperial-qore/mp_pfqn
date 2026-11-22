#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "pfqn_sing.h"

/* Shared singularity oracle for the MoM / CoMoM / gmom family.
 *
 * The per-level coefficient matrix of every member is built from the
 * loadings of the NON-recursion classes 1..R-1 only; the recursion class
 * R enters the right-hand side, never the coefficient matrix (paper
 * Sec 5.1: the coefficient matrix is independent of the population being
 * processed, and of the recursion class's demand).  Hence:
 *
 *   1. singularity is a property of the demands alone, decidable ONCE at
 *      initialisation without running the recursion; and
 *   2. if a single class causes the degeneracy, MOVING IT TO THE
 *      RECURSION POSITION removes it, because its loadings then no longer
 *      appear in the coefficient matrix.
 *
 * The oracle builds the b=1 generalized-MoM coefficient matrix A via
 * setup1 (whose entries depend only on L, not on N or Z), forms the exact
 * normal-equations matrix A^T A, and tests its rank by fraction-free LU.
 * This is the shared demand-degeneracy indicator for the whole family:
 * the members differ in basis but share the loading relations that make
 * the CE/PC system rank-deficient.
 */

/* Test whether the recursion coefficient matrix is singular for the model
 * with demand matrix L (M queues x R classes), with class R-1 (the last
 * column) as the recursion class.  Returns 1 if singular, 0 if full rank.
 * Population/think-time independent, so N and Z are passed as zeros. */
static int level_singular(mpz_t** L, int m, int r)
{
	int i, s;
	int* N = (int*) calloc(r, sizeof(int));
	mpz_t* Z = (mpz_t*) malloc(r * sizeof(mpz_t));
	for (s = 0; s < r; s++) mpz_init(Z[s]);

	LS1* ls = setup1(L, m, r, N, Z);
	int colA = ls->colA;
	mpq_mat_t AtA = mpq_matzeros(colA, colA);
	mpq_mattransmul(AtA, ls->A, ls->A, ls->nrows, colA);
	int* idx = mpq_ludcmp(AtA, colA);
	int singular = (idx == NULL);
	if (idx) free(idx);
	for (i = 0; i < colA+1; i++) { for (s = 0; s < colA+1; s++) mpq_clear(AtA[i][s]); free(AtA[i]); }
	free(AtA);
	ls1_free(ls, r);
	for (s = 0; s < r; s++) mpz_clear(Z[s]);
	free(Z); free(N);
	return singular;
}

int pfqn_recursion_singular(mpz_t** L, int M, int R)
{
	/* The recursion builds a coefficient matrix at every level (r, m),
	 * 2 <= r <= R classes processed, 2 <= m <= M queues in the prefix
	 * sub-model.  All are loading-only (population/think-time
	 * independent).  The recursion is degenerate iff ANY of them is
	 * singular, so test them all. */
	int r, m;
	if (R < 2) return 0;
	for (r = 2; r <= R; r++)
		for (m = 2; m <= M; m++)
			if (level_singular(L, m, r)) return 1;
	return 0;
}

/* Build a column-permuted copy of L placing original class `last` in the
 * final (recursion) column, keeping the other classes in order. */
static mpz_t** perm_last(mpz_t** L, int M, int R, int last)
{
	int k, s, c;
	mpz_t** P = (mpz_t**) malloc(M * sizeof(mpz_t*));
	for (k = 0; k < M; k++) {
		P[k] = (mpz_t*) malloc(R * sizeof(mpz_t));
		c = 0;
		for (s = 0; s < R; s++) { if (s == last) continue; mpz_init_set(P[k][c++], L[k][s]); }
		mpz_init_set(P[k][R-1], L[k][last]);
	}
	return P;
}
static void perm_free(mpz_t** P, int M, int R)
{
	int k, s;
	for (k = 0; k < M; k++) { for (s = 0; s < R; s++) mpz_clear(P[k][s]); free(P[k]); }
	free(P);
}

/* Find a class that, used as the recursion class (moved to the last
 * column), makes the coefficient matrix non-singular.  Returns that
 * original class index in 0..R-1, or -1 if no single-class reorder
 * removes the degeneracy (a genuine loading relation needing the
 * convolution fallback).  Tries the current last class first. */
int pfqn_nonsingular_recclass(mpz_t** L, int M, int R)
{
	int c;
	if (!pfqn_recursion_singular(L, M, R)) return R-1;   /* already fine */
	for (c = 0; c < R; c++) {
		if (c == R-1) continue;
		mpz_t** P = perm_last(L, M, R, c);
		int sing = pfqn_recursion_singular(P, M, R);
		perm_free(P, M, R);
		if (!sing) return c;
	}
	return -1;
}
