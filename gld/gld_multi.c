#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gld.h"

/*
 * Multi-class load-dependent normalizing constant via recursion.
 * Port of pfqn_gld.m from LINE-dev.
 *
 * Base cases:
 *   Nt == 0 -> G = 1
 *   M == 0  -> G = 0
 *   R == 1  -> call gld_single
 *   M == 1  -> multinomial formula
 * Recursion:
 *   G = gld(L[0:M-2], N, mu[0:M-2], M-1, R, Nt)
 *     + sum_r [L[M-1][r]/mu[M-1][0]] * gld(L, N-e_r, mushift(mu,M-1), M, R, Nt-1)
 */

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */
void gld_auto_mu(qnmodel *qn)
{
	int m, k;
	if (qn->mu != NULL || qn->Nt <= 0)
		return;
	qn->mu = (mpq_t **)malloc(qn->M * sizeof(mpq_t *));
	for (m = 0; m < qn->M; m++) {
		qn->mu[m] = (mpq_t *)malloc(qn->Nt * sizeof(mpq_t));
		for (k = 0; k < qn->Nt; k++) {
			mpq_init(qn->mu[m][k]);
			mpq_set_ui(qn->mu[m][k], MIN(k + 1, qn->mi[m]), 1);
		}
	}
	qn->isLD = 1;
}

void gld_multi(mpq_t G, mpz_t **L, int *N, mpq_t **mu, int M, int R, int Nt)
{
	int m, r, k;

	/* Base case: all populations zero */
	if (Nt == 0) {
		mpq_set_ui(G, 1, 1);
		return;
	}

	/* Base case: no stations */
	if (M == 0) {
		mpq_set_ui(G, 0, 1);
		return;
	}

	/* Base case: single class -> use efficient DP */
	if (R == 1) {
		mpz_t *L_single = (mpz_t *)malloc(M * sizeof(mpz_t));
		for (m = 0; m < M; m++)
			mpz_init_set(L_single[m], L[m][0]);
		gld_single(G, L_single, N[0], mu, M);
		for (m = 0; m < M; m++)
			mpz_clear(L_single[m]);
		free(L_single);
		return;
	}

	/* Base case: single station -> multinomial formula
	 * G = Nt! / prod(N[r]!) * prod(L[0][r]^N[r]) / prod(mu[0][k])
	 */
	if (M == 1) {
		/* Check if any class has positive pop but zero demand */
		for (r = 0; r < R; r++) {
			if (N[r] > 0 && mpz_sgn(L[0][r]) == 0) {
				mpq_set_ui(G, 0, 1);
				return;
			}
		}

		/* Numerator: Nt! * prod(L[0][r]^N[r]) */
		mpz_t num, L_power, factval;
		mpz_init(num);
		mpz_init(L_power);
		mpz_init(factval);

		mpz_fac_ui(num, Nt);
		for (r = 0; r < R; r++) {
			if (N[r] > 0 && mpz_sgn(L[0][r]) > 0) {
				mpz_pow_ui(L_power, L[0][r], N[r]);
				mpz_mul(num, num, L_power);
			}
		}

		/* Denominator: prod(N[r]!) */
		mpz_t den_int;
		mpz_init(den_int);
		mpz_set_ui(den_int, 1);
		for (r = 0; r < R; r++) {
			mpz_fac_ui(factval, N[r]);
			mpz_mul(den_int, den_int, factval);
		}

		/* prod(mu[0][k]) for k=0..Nt-1 */
		mpq_t mu_prod;
		mpq_init(mu_prod);
		mpq_set_ui(mu_prod, 1, 1);
		for (k = 0; k < Nt; k++)
			mpq_mul(mu_prod, mu_prod, mu[0][k]);

		/* G = num / (den_int * mu_prod) */
		mpq_t result, den_q;
		mpq_init(result);
		mpq_init(den_q);
		mpq_set_z(result, num);
		mpq_set_z(den_q, den_int);
		mpq_mul(den_q, den_q, mu_prod);
		mpq_div(result, result, den_q);
		mpq_set(G, result);

		mpz_clear(num);
		mpz_clear(L_power);
		mpz_clear(factval);
		mpz_clear(den_int);
		mpq_clear(mu_prod);
		mpq_clear(result);
		mpq_clear(den_q);
		return;
	}

	/* General recursion */
	mpq_t G_sub, term, ratio;
	mpq_init(G_sub);
	mpq_init(term);
	mpq_init(ratio);

	/* First term: strip last station -> gld(L, N, mu, M-1, R, Nt) */
	gld_multi(G_sub, L, N, mu, M - 1, R, Nt);
	mpq_set(G, G_sub);

	/* Second term: for each class r with N[r] > 0 */
	for (r = 0; r < R; r++) {
		if (N[r] > 0) {
			/* ratio = L[M-1][r] / mu[M-1][0] */
			mpq_set_z(ratio, L[M - 1][r]);
			mpq_div(ratio, ratio, mu[M - 1][0]);

			/* N_sub = N with N[r] decremented */
			int *N_sub = (int *)malloc(R * sizeof(int));
			for (int s = 0; s < R; s++)
				N_sub[s] = N[s];
			N_sub[r]--;

			if (Nt - 1 == 0) {
				/* All populations zero -> G_sub = 1 */
				mpq_set_ui(G_sub, 1, 1);
			} else {
				/* mu_shifted = mushift(mu, M-1): shift station M-1 left, keep others */
				mpq_t **mu_shifted = (mpq_t **)malloc(M * sizeof(mpq_t *));
				for (m = 0; m < M; m++) {
					mu_shifted[m] = (mpq_t *)malloc((Nt - 1) * sizeof(mpq_t));
					for (k = 0; k < Nt - 1; k++) {
						mpq_init(mu_shifted[m][k]);
						if (m == M - 1)
							mpq_set(mu_shifted[m][k], mu[m][k + 1]);
						else
							mpq_set(mu_shifted[m][k], mu[m][k]);
					}
				}

				gld_multi(G_sub, L, N_sub, mu_shifted, M, R, Nt - 1);

				for (m = 0; m < M; m++) {
					for (k = 0; k < Nt - 1; k++)
						mpq_clear(mu_shifted[m][k]);
					free(mu_shifted[m]);
				}
				free(mu_shifted);
			}

			/* G += ratio * G_sub */
			mpq_mul(term, ratio, G_sub);
			mpq_add(G, G, term);

			free(N_sub);
		}
	}

	mpq_clear(G_sub);
	mpq_clear(term);
	mpq_clear(ratio);
}
