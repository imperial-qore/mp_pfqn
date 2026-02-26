#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gld.h"

/*
 * Single-class load-dependent normalizing constant via 3D DP.
 * Port of pfqn_gldsingle.m from LINE-dev.
 *
 * g(m, n, tm) for m=0..M, n=0..N, tm=0..N+1
 * Base cases:
 *   g(0, n, 1) = 0 for n>=1
 *   g(m, 0, tm) = 1 for m>=1, tm>=1
 * Recursion:
 *   g(m, n, tm) = g(m-1, n, 1) + L[m-1]*g(m, n-1, tm+1)/mu[m-1][tm-1]
 * Result: G = g(M, N, 1)
 */
void gld_single(mpq_t G, mpz_t *L, int N, mpq_t **mu, int M)
{
	int m, n, tm;

	if (N == 0) {
		mpq_set_ui(G, 1, 1);
		return;
	}
	if (M == 0) {
		mpq_set_ui(G, 0, 1);
		return;
	}

	/* Allocate 3D array g[M+1][N+1][N+2] */
	mpq_t ***g = (mpq_t ***)malloc((M + 1) * sizeof(mpq_t **));
	for (m = 0; m <= M; m++) {
		g[m] = (mpq_t **)malloc((N + 1) * sizeof(mpq_t *));
		for (n = 0; n <= N; n++) {
			g[m][n] = (mpq_t *)malloc((N + 2) * sizeof(mpq_t));
			for (tm = 0; tm < N + 2; tm++)
				mpq_init(g[m][n][tm]);
		}
	}

	/* Base case: g[0][n][1] = 0 for n=1..N (already 0 from init) */
	/* Base case: g[m][0][tm] = 1 for m=1..M, tm=1..N+1 */
	for (m = 1; m <= M; m++)
		for (tm = 1; tm <= N + 1; tm++)
			mpq_set_ui(g[m][0][tm], 1, 1);

	/* Recursion */
	mpq_t term;
	mpq_init(term);
	for (m = 1; m <= M; m++) {
		for (n = 1; n <= N; n++) {
			for (tm = 1; tm <= N - n + 1; tm++) {
				/* g[m][n][tm] = g[m-1][n][1] + L[m-1]*g[m][n-1][tm+1]/mu[m-1][tm-1] */
				mpq_set(g[m][n][tm], g[m - 1][n][1]);
				mpq_set_z(term, L[m - 1]);
				mpq_mul(term, term, g[m][n - 1][tm + 1]);
				mpq_div(term, term, mu[m - 1][tm - 1]);
				mpq_add(g[m][n][tm], g[m][n][tm], term);
			}
		}
	}
	mpq_clear(term);

	mpq_set(G, g[M][N][1]);

	/* Free */
	for (m = 0; m <= M; m++) {
		for (n = 0; n <= N; n++) {
			for (tm = 0; tm < N + 2; tm++)
				mpq_clear(g[m][n][tm]);
			free(g[m][n]);
		}
		free(g[m]);
	}
	free(g);
}
