#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "comomld.h"

/*
 * CoMoM-LD: Class-Oriented Method of Moments for repairman models
 * with load-dependent service rates.
 *
 * Port of pfqn_comomrm_ld.m from LINE-dev.
 *
 * Algorithm:
 * 1. Detect infinite servers: station m is IS if mu[m][k] == k+1 for all k.
 *    Absorb IS demands into Z (think times).
 * 2. After absorption, exactly M_q == 1 queueing station must remain.
 * 3. Class-by-class bidiagonal recursion:
 *    - Tr is bidiagonal (Nt+1)x(Nt+1):
 *        diagonal = Z_eff[r]
 *        superdiagonal[i] = L_eff[r] * (Nt - i) / mu_q[Nt - i - 1]  for i=0..Nt-1
 *    - h = [0, ..., 0, 1] (length Nt+1, last element = 1)
 *    - For each class r, for nr=1..N[r]: h = (1/nr) * Tr * h, then
 *      scale = sum(h), h = h / scale
 * 4. G = product of all scale values
 * 5. prob[k] = h[Nt - k] for k=0..Nt
 */
void comomrm_ld(mpq_t G, mpq_t *prob, qnmodel *qn)
{
	int M = qn->M;
	int R = qn->R;
	int Nt = qn->Nt;
	int m, r, k, i;

	/* --- Step 1: Detect infinite servers and absorb into Z --- */

	/* is_inf[m] = 1 if station m is an infinite server */
	int *is_inf = (int *)calloc(M, sizeof(int));
	int num_inf = 0;

	for (m = 0; m < M; m++) {
		int all_match = 1;
		for (k = 0; k < Nt; k++) {
			/* Check if mu[m][k] == k+1 (as rational) */
			mpq_t expected;
			mpq_init(expected);
			mpq_set_ui(expected, k + 1, 1);
			if (!mpq_equal(qn->mu[m][k], expected)) {
				all_match = 0;
				mpq_clear(expected);
				break;
			}
			mpq_clear(expected);
		}
		if (all_match) {
			is_inf[m] = 1;
			num_inf++;
		}
	}

	int M_q = M - num_inf;

	/* Compute effective Z: Z_eff[r] = Z[r] + sum of L[m][r] over IS stations */
	mpq_t *Z_eff = (mpq_t *)malloc(R * sizeof(mpq_t));
	for (r = 0; r < R; r++) {
		mpq_init(Z_eff[r]);
		mpq_set_z(Z_eff[r], qn->Z[r]);
		for (m = 0; m < M; m++) {
			if (is_inf[m]) {
				mpq_t Lval;
				mpq_init(Lval);
				mpq_set_z(Lval, qn->L[m][r]);
				mpq_add(Z_eff[r], Z_eff[r], Lval);
				mpq_clear(Lval);
			}
		}
	}

	/* Find the single queueing station index */
	int q_idx = -1;
	if (M_q == 0) {
		/* No queueing stations: model has only delays.
		 * G = prod_r Z_eff[r]^N[r] / prod_r N[r]!  (multinomial of delays)
		 * Actually for pure delay network:
		 * G = Nt! / prod(N[r]!) * prod(Z_eff[r]^N[r])
		 * All jobs are always at the delay station, prob[Nt] = 1 (if station existed)
		 * But there's no queueing station, so prob[0] = 1.
		 */
		mpz_t num, factval;
		mpz_init(num);
		mpz_init(factval);

		/* G_num = Nt! * prod(Z_eff[r]^N[r]) -- but Z_eff are rationals */
		mpq_t G_val, Z_pow;
		mpq_init(G_val);
		mpq_init(Z_pow);

		/* Nt! */
		mpz_fac_ui(num, Nt);
		mpq_set_z(G_val, num);

		for (r = 0; r < R; r++) {
			/* Z_eff[r]^N[r] */
			mpq_set_ui(Z_pow, 1, 1);
			for (int j = 0; j < qn->N[r]; j++)
				mpq_mul(Z_pow, Z_pow, Z_eff[r]);
			mpq_mul(G_val, G_val, Z_pow);
		}

		/* divide by prod(N[r]!) */
		mpz_t den;
		mpz_init(den);
		mpz_set_ui(den, 1);
		for (r = 0; r < R; r++) {
			mpz_fac_ui(factval, qn->N[r]);
			mpz_mul(den, den, factval);
		}
		mpq_t den_q;
		mpq_init(den_q);
		mpq_set_z(den_q, den);
		mpq_div(G_val, G_val, den_q);

		mpq_set(G, G_val);

		/* prob: all jobs at delay, none at queueing station */
		for (k = 0; k <= Nt; k++)
			mpq_set_ui(prob[k], 0, 1);
		mpq_set_ui(prob[0], 1, 1);

		mpz_clear(num);
		mpz_clear(factval);
		mpz_clear(den);
		mpq_clear(G_val);
		mpq_clear(Z_pow);
		mpq_clear(den_q);
		for (r = 0; r < R; r++)
			mpq_clear(Z_eff[r]);
		free(Z_eff);
		free(is_inf);
		return;
	}

	if (M_q != 1) {
		fprintf(stderr, "comomld: requires exactly 1 queueing station after IS absorption, found %d\n", M_q);
		mpq_set_ui(G, 0, 1);
		for (k = 0; k <= Nt; k++)
			mpq_set_ui(prob[k], 0, 1);
		for (r = 0; r < R; r++)
			mpq_clear(Z_eff[r]);
		free(Z_eff);
		free(is_inf);
		return;
	}

	for (m = 0; m < M; m++) {
		if (!is_inf[m]) {
			q_idx = m;
			break;
		}
	}

	/* Effective L and mu for the single queueing station */
	/* L_eff[r] = L[q_idx][r], mu_eff[k] = mu[q_idx][k] */

	/* --- Step 3: Bidiagonal recursion --- */

	/* h vector of size Nt+1: h[Nt] = 1, rest = 0 */
	mpq_t *h = (mpq_t *)malloc((Nt + 1) * sizeof(mpq_t));
	mpq_t *h_new = (mpq_t *)malloc((Nt + 1) * sizeof(mpq_t));
	for (i = 0; i <= Nt; i++) {
		mpq_init(h[i]);
		mpq_init(h_new[i]);
	}
	mpq_set_ui(h[Nt], 1, 1);

	/* Accumulate G = product of scale factors */
	mpq_set_ui(G, 1, 1);

	mpq_t term, scale, inv_nr, L_q_r;
	mpq_init(term);
	mpq_init(scale);
	mpq_init(inv_nr);
	mpq_init(L_q_r);

	for (r = 0; r < R; r++) {
		/* L_q_r = L[q_idx][r] as rational */
		mpq_set_z(L_q_r, qn->L[q_idx][r]);

		/* Build Tr and apply: h_new = (1/nr) * Tr * h
		 * Tr is bidiagonal:
		 *   Tr[i][i] = Z_eff[r]              for i = 0..Nt
		 *   Tr[i][i+1] = L_q_r * (Nt - i) / mu_q[Nt - i - 1]  for i = 0..Nt-1
		 *
		 * So (Tr * h)[i] = Z_eff[r] * h[i] + superdiag[i] * h[i+1]
		 *   for i < Nt
		 * (Tr * h)[Nt] = Z_eff[r] * h[Nt]
		 *
		 * This is O(Nt) per step.
		 */
		for (int nr = 1; nr <= qn->N[r]; nr++) {
			/* inv_nr = 1/nr */
			mpq_set_ui(inv_nr, 1, nr);

			/* Compute h_new = (1/nr) * Tr * h */

			/* Last element: h_new[Nt] = (1/nr) * Z_eff[r] * h[Nt] */
			mpq_mul(h_new[Nt], Z_eff[r], h[Nt]);
			mpq_mul(h_new[Nt], inv_nr, h_new[Nt]);

			/* Elements i = 0..Nt-1:
			 * h_new[i] = (1/nr) * (Z_eff[r] * h[i] + L_q_r * (Nt-i) / mu_q[Nt-i-1] * h[i+1])
			 */
			for (i = Nt - 1; i >= 0; i--) {
				/* diagonal part: Z_eff[r] * h[i] */
				mpq_mul(h_new[i], Z_eff[r], h[i]);

				/* superdiagonal part: L_q_r * (Nt - i) / mu_q[Nt-i-1] * h[i+1] */
				mpq_set_z(term, qn->L[q_idx][r]); /* L_q_r as fresh rational */
				mpq_mul(term, term, h[i + 1]);

				/* multiply by (Nt - i) */
				mpq_t factor;
				mpq_init(factor);
				mpq_set_ui(factor, Nt - i, 1);
				mpq_mul(term, term, factor);

				/* divide by mu_q[Nt - i - 1] */
				mpq_div(term, term, qn->mu[q_idx][Nt - i - 1]);

				mpq_add(h_new[i], h_new[i], term);

				/* multiply by 1/nr */
				mpq_mul(h_new[i], inv_nr, h_new[i]);

				mpq_clear(factor);
			}

			/* Compute scale = sum(h_new) */
			mpq_set_ui(scale, 0, 1);
			for (i = 0; i <= Nt; i++)
				mpq_add(scale, scale, h_new[i]);

			/* G *= scale */
			mpq_mul(G, G, scale);

			/* h = h_new / scale */
			for (i = 0; i <= Nt; i++) {
				mpq_div(h[i], h_new[i], scale);
			}
		}
	}

	/* --- Step 5: prob[k] = h[Nt - k] --- */
	/* h is already normalized (sum = 1) */
	for (k = 0; k <= Nt; k++)
		mpq_set(prob[k], h[Nt - k]);

	/* Cleanup */
	mpq_clear(term);
	mpq_clear(scale);
	mpq_clear(inv_nr);
	mpq_clear(L_q_r);

	for (i = 0; i <= Nt; i++) {
		mpq_clear(h[i]);
		mpq_clear(h_new[i]);
	}
	free(h);
	free(h_new);

	for (r = 0; r < R; r++)
		mpq_clear(Z_eff[r]);
	free(Z_eff);
	free(is_inf);
}
