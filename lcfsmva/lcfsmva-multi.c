#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

double lcfsmva_multi(qnmodel* qn, mpq_t *X, mpq_t **Q, mpq_t G, mpq_t **B_out)
{
	struct rusage ruse;
	double t_start = CPUTIME;
	int R = qn->R;
	int r, k, s;

	/* alpha[r] = L[0][r], beta[r] = L[1][r] — service demands as mpz_t */
	/* These are just aliases, no need to copy */
	mpz_t *alpha = qn->L[0]; /* station 0: LCFS-NP */
	mpz_t *beta  = qn->L[1]; /* station 1: LCFS-PR */

	int *n = (int*)initpop(R);
	if (n == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory in lcfsmva_multi\n");
		return -1.0;
	}
	int *planesizes = getplanesizes(qn->N, R);
	if (planesizes == NULL) {
		free(n);
		fprintf(stderr, "Error: Failed to allocate memory in lcfsmva_multi\n");
		return -1.0;
	}
	int nstates = planesizes[0];

	/* Allocate B[nstates][2][R], QN[nstates][2][R], TN[nstates][R] */
	mpq_t ***B = (mpq_t***)calloc(nstates, sizeof(mpq_t**));
	mpq_t ***QN = (mpq_t***)calloc(nstates, sizeof(mpq_t**));
	mpq_t **TN = (mpq_t**)calloc(nstates, sizeof(mpq_t*));
	if (B == NULL || QN == NULL || TN == NULL) {
		fprintf(stderr, "Error: Failed to allocate state arrays (%d states). Model may be too large.\n", nstates);
		free(n); free(planesizes);
		if (B) free(B);
		if (QN) free(QN);
		if (TN) free(TN);
		return -1.0;
	}

	for (int i = 0; i < nstates; i++) {
		B[i] = (mpq_t**)calloc(2, sizeof(mpq_t*));
		QN[i] = (mpq_t**)calloc(2, sizeof(mpq_t*));
		TN[i] = (mpq_t*)calloc(R, sizeof(mpq_t));
		if (B[i] == NULL || QN[i] == NULL || TN[i] == NULL) {
			fprintf(stderr, "Error: Failed to allocate state %d. Model may be too large.\n", i);
			/* Partial cleanup */
			for (int j = 0; j <= i; j++) {
				if (TN[j]) {
					for (r = 0; r < R; r++) mpq_clear(TN[j][r]);
					free(TN[j]);
				}
				if (QN[j]) {
					for (s = 0; s < 2; s++) {
						if (QN[j][s]) {
							for (r = 0; r < R; r++) mpq_clear(QN[j][s][r]);
							free(QN[j][s]);
						}
					}
					free(QN[j]);
				}
				if (B[j]) {
					for (s = 0; s < 2; s++) {
						if (B[j][s]) {
							for (r = 0; r < R; r++) mpq_clear(B[j][s][r]);
							free(B[j][s]);
						}
					}
					free(B[j]);
				}
			}
			free(B); free(QN); free(TN);
			free(n); free(planesizes);
			return -1.0;
		}
		for (s = 0; s < 2; s++) {
			B[i][s] = (mpq_t*)calloc(R, sizeof(mpq_t));
			QN[i][s] = (mpq_t*)calloc(R, sizeof(mpq_t));
			if (B[i][s] == NULL || QN[i][s] == NULL) {
				fprintf(stderr, "Error: Failed to allocate station arrays for state %d.\n", i);
				/* Full cleanup would be complex; exit */
				free(n); free(planesizes);
				return -1.0;
			}
			for (r = 0; r < R; r++) {
				mpq_init(B[i][s][r]);
				mpq_init(QN[i][s][r]);
			}
		}
		for (r = 0; r < R; r++)
			mpq_init(TN[i][r]);
	}

	/* Allocate CA arrays for normalizing constant: G_ca[nstates], V[nstates] */
	mpq_t *G_ca = (mpq_t*)calloc(nstates, sizeof(mpq_t));
	mpq_t *V = (mpq_t*)calloc(nstates, sizeof(mpq_t));
	if (G_ca == NULL || V == NULL) {
		fprintf(stderr, "Error: Failed to allocate G_ca/V arrays.\n");
		free(n); free(planesizes);
		if (G_ca) free(G_ca);
		if (V) free(V);
		return -1.0;
	}
	for (int i = 0; i < nstates; i++) {
		mpq_init(G_ca[i]);
		mpq_init(V[i]);
	}

	/* Temporary GMP variables used in the main loop */
	mpq_t tmp_q, tmp_q2, tmp_q3, tmp_q4;
	mpz_t tmp_z, tmp_z2;
	mpq_init(tmp_q);
	mpq_init(tmp_q2);
	mpq_init(tmp_q3);
	mpq_init(tmp_q4);
	mpz_init(tmp_z);
	mpz_init(tmp_z2);

	/* scale factors */
	mpz_t scale_np, scale_pr_k;
	mpz_init(scale_np);
	mpz_init(scale_pr_k);

	/* For the CA denom computation */
	mpq_t denom_q;
	mpq_init(denom_q);

	int curindex = -1;

	/* Population iteration */
	do
	{
		curindex++;

		/* Progress output */
		int* widths = (int*)calloc(R, sizeof(int));
		for (r = 0; r < R; r++) {
			int val = qn->N[r];
			widths[r] = 1;
			while (val >= 10) { widths[r]++; val /= 10; }
		}
		fprintf(stderr, "\rn=(");
		for (r = 0; r < R; r++) {
			fprintf(stderr, "%*d", widths[r], n[r]);
			if (r < R - 1) fprintf(stderr, ",");
		}
		double t_current = CPUTIME;
		fprintf(stderr, ") - Time: %.2f s  ", t_current - t_start);
		fflush(stderr);
		free(widths);

		int idx = popindex(n, R, planesizes);
		int sumn = sum(n, R);

		if (sumn == 0) {
			/* Base: G_ca[0] = 1, V[0] = 1 */
			mpq_set_ui(G_ca[idx], 1, 1);
			mpq_set_ui(V[idx], 1, 1);
			/* B, QN, TN already zeroed by mpq_init */
			continue;
		}

		/* ============ CA normalizing constant recursion ============ */
		/* V[n] = prod(alpha[r]^n[r]) * sum_{r: n[r]>0} V[n-e_r] */
		/* G_ca[n] = V[n] + sum_{r: n[r]>0} alpha[r]^(sumn-1) * beta[r] * G_ca[n-e_r] */

		/* Compute prod(alpha[r]^n[r]) into tmp_z */
		mpz_set_ui(tmp_z, 1);
		for (r = 0; r < R; r++) {
			if (n[r] > 0) {
				mpz_pow_ui(tmp_z2, alpha[r], (unsigned long)n[r]);
				mpz_mul(tmp_z, tmp_z, tmp_z2);
			}
		}
		/* tmp_z now holds prod(alpha[r]^n[r]) */

		/* sum_{r: n[r]>0} V[n-e_r] */
		mpq_set_ui(tmp_q, 0, 1);
		for (r = 0; r < R; r++) {
			if (n[r] > 0) {
				n[r]--;
				int idx_r = popindex(n, R, planesizes);
				n[r]++;
				mpq_add(tmp_q, tmp_q, V[idx_r]);
			}
		}
		/* V[n] = tmp_z * tmp_q */
		mpq_set_z(tmp_q2, tmp_z);
		mpq_mul(V[idx], tmp_q2, tmp_q);

		/* G_ca[n] = V[n] + sum_{r: n[r]>0} alpha[r]^(sumn-1) * beta[r] * G_ca[n-e_r] */
		mpq_set(G_ca[idx], V[idx]);
		for (r = 0; r < R; r++) {
			if (n[r] > 0) {
				n[r]--;
				int idx_r = popindex(n, R, planesizes);
				n[r]++;
				/* alpha[r]^(sumn-1) * beta[r] */
				mpz_pow_ui(tmp_z2, alpha[r], (unsigned long)(sumn - 1));
				mpz_mul(tmp_z2, tmp_z2, beta[r]);
				mpq_set_z(tmp_q, tmp_z2);
				mpq_mul(tmp_q, tmp_q, G_ca[idx_r]);
				mpq_add(G_ca[idx], G_ca[idx], tmp_q);
			}
		}

		/* ============ LCFS MVA recursion ============ */

		if (sumn == 1) {
			/* Base case: single job of class r (where n[r]==1) */
			for (r = 0; r < R; r++) {
				if (n[r] == 1) {
					/* denom = alpha[r] + beta[r] */
					mpz_add(tmp_z, alpha[r], beta[r]);
					/* Q[idx][0][r] = alpha[r] / denom */
					mpq_set_z(QN[idx][0][r], alpha[r]);
					mpq_set_z(tmp_q, tmp_z);
					mpq_div(QN[idx][0][r], QN[idx][0][r], tmp_q);
					/* Q[idx][1][r] = beta[r] / denom */
					mpq_set_z(QN[idx][1][r], beta[r]);
					mpq_div(QN[idx][1][r], QN[idx][1][r], tmp_q);
					/* B[idx][0][r] = alpha[r] / denom */
					mpq_set(B[idx][0][r], QN[idx][0][r]);
					/* B[idx][1][r] = beta[r] / denom */
					mpq_set(B[idx][1][r], QN[idx][1][r]);
					/* T[idx][r] = 1 / denom */
					mpq_set_ui(TN[idx][r], 1, 1);
					mpq_div(TN[idx][r], TN[idx][r], tmp_q);
					break; /* Only one class has n[r]==1 */
				}
			}
		} else {
			/* Recursive case: sumn > 1 */

			/* 1. Compute scale_np = prod_r alpha[r]^n[r] */
			mpz_set_ui(scale_np, 1);
			for (r = 0; r < R; r++) {
				if (n[r] > 0) {
					mpz_pow_ui(tmp_z, alpha[r], (unsigned long)n[r]);
					mpz_mul(scale_np, scale_np, tmp_z);
				}
			}

			/* For each class k with n[k] > 0 */
			for (k = 0; k < R; k++) {
				if (n[k] == 0) continue;

				/* scale_pr_k = alpha[k]^(sumn-1) * beta[k] */
				mpz_pow_ui(scale_pr_k, alpha[k], (unsigned long)(sumn - 1));
				mpz_mul(scale_pr_k, scale_pr_k, beta[k]);

				/* 2. Compute unscaled waiting times */
				/* Wnp = 1 + QN[n-e_k][0][k] */
				n[k]--;
				int idx_k = popindex(n, R, planesizes);
				n[k]++;

				mpq_t Wnp, Wpr;
				mpq_init(Wnp);
				mpq_init(Wpr);
				mpq_set_ui(Wnp, 1, 1);
				mpq_add(Wnp, Wnp, QN[idx_k][0][k]);
				/* Wpr = 1 + QN[n-e_k][1][k] */
				mpq_set_ui(Wpr, 1, 1);
				mpq_add(Wpr, Wpr, QN[idx_k][1][k]);

				/* For each class r != k with n[r] > 0 */
				for (r = 0; r < R; r++) {
					if (r == k || n[r] == 0) continue;

					n[r]--;
					int idx_r = popindex(n, R, planesizes);
					n[r]++;

					/* LCFS-NP contribution */
					/* if B[n-e_r][0][k] != 0: */
					if (mpq_sgn(B[idx_r][0][k]) != 0) {
						/* ratio_np = B[n-e_k][0][r] / B[n-e_r][0][k] */
						mpq_div(tmp_q, B[idx_k][0][r], B[idx_r][0][k]);
						/* (alpha[k]/alpha[r]) * ratio_np * QN[n-e_r][0][k] */
						mpq_set_z(tmp_q2, alpha[k]);
						mpq_set_z(tmp_q3, alpha[r]);
						mpq_div(tmp_q2, tmp_q2, tmp_q3);
						mpq_mul(tmp_q2, tmp_q2, tmp_q);
						mpq_mul(tmp_q2, tmp_q2, QN[idx_r][0][k]);
						mpq_add(Wnp, Wnp, tmp_q2);
					}

					/* LCFS-PR contribution */
					/* if B[n-e_r][1][k] != 0: */
					if (mpq_sgn(B[idx_r][1][k]) != 0) {
						/* ratio_pr = B[n-e_k][1][r] / B[n-e_r][1][k] */
						mpq_div(tmp_q, B[idx_k][1][r], B[idx_r][1][k]);
						/* (alpha[r]/alpha[k]) * ratio_pr * QN[n-e_r][1][k] */
						mpq_set_z(tmp_q2, alpha[r]);
						mpq_set_z(tmp_q3, alpha[k]);
						mpq_div(tmp_q2, tmp_q2, tmp_q3);
						mpq_mul(tmp_q2, tmp_q2, tmp_q);
						mpq_mul(tmp_q2, tmp_q2, QN[idx_r][1][k]);
						mpq_add(Wpr, Wpr, tmp_q2);
					}
				}

				/* 3. Back probabilities */
				/* Wnp_total = scale_np * Wnp, Wpr_total = scale_pr_k * Wpr */
				/* denom = scale_np * Wnp + scale_pr_k * Wpr */
				mpq_set_z(tmp_q, scale_np);
				mpq_mul(tmp_q, tmp_q, Wnp);   /* tmp_q = scale_np * Wnp */

				mpq_set_z(tmp_q2, scale_pr_k);
				mpq_mul(tmp_q2, tmp_q2, Wpr);  /* tmp_q2 = scale_pr_k * Wpr */

				mpq_add(denom_q, tmp_q, tmp_q2); /* denom = tmp_q + tmp_q2 */

				/* B[idx][0][k] = n[k] * scale_np / denom */
				/* But we need n[k] * scale_np * Wnp / denom? No.
				   Actually: B[idx][0][k] = n[k] * scale_np / denom
				   where denom = scale_np * Wnp + scale_pr_k * Wpr
				   Wait, the spec says:
				   B[idx][0][k] = n[k] * scale_np / denom
				   This uses the raw scale_np (not multiplied by Wnp).
				   Let me re-read...

				   "denom = scale_np * Wnp + scale_pr_k * Wpr"
				   "B[idx][0][k] = n[k] * scale_np / denom"

				   So denom already includes Wnp/Wpr, and B uses raw scale_np.
				*/
				/* Recompute: denom = scale_np * Wnp + scale_pr_k * Wpr */
				/* Already in denom_q */

				/* B[idx][0][k] = n[k] * scale_np / denom */
				mpq_set_z(tmp_q3, scale_np);
				mpq_set_ui(tmp_q4, (unsigned long)n[k], 1);
				mpq_mul(tmp_q3, tmp_q3, tmp_q4);
				mpq_div(B[idx][0][k], tmp_q3, denom_q);

				/* B[idx][1][k] = n[k] * scale_pr_k / denom */
				mpq_set_z(tmp_q3, scale_pr_k);
				mpq_mul(tmp_q3, tmp_q3, tmp_q4);
				mpq_div(B[idx][1][k], tmp_q3, denom_q);

				mpq_clear(Wnp);
				mpq_clear(Wpr);
			}

			/* 4. Queue lengths */
			for (k = 0; k < R; k++) {
				if (n[k] == 0) continue;
				for (s = 0; s < 2; s++) {
					/* QN[idx][s][k] = B[idx][s][k] */
					mpq_set(QN[idx][s][k], B[idx][s][k]);
					/* + sum_r B[idx][s][r] * QN[n-e_r][s][k] */
					for (r = 0; r < R; r++) {
						if (n[r] == 0) continue;
						n[r]--;
						int idx_r = popindex(n, R, planesizes);
						n[r]++;
						mpq_mul(tmp_q, B[idx][s][r], QN[idx_r][s][k]);
						mpq_add(QN[idx][s][k], QN[idx][s][k], tmp_q);
					}
				}
			}

			/* 5. Throughputs */
			for (k = 0; k < R; k++) {
				if (n[k] == 0) continue;

				n[k]--;
				int idx_k = popindex(n, R, planesizes);
				n[k]++;

				/* Unp = sum_r alpha[r] * TN[idx_k][r] */
				mpq_t Unp;
				mpq_init(Unp);
				mpq_set_ui(Unp, 0, 1);
				for (r = 0; r < R; r++) {
					mpq_set_z(tmp_q, alpha[r]);
					mpq_mul(tmp_q, tmp_q, TN[idx_k][r]);
					mpq_add(Unp, Unp, tmp_q);
				}

				/* TN[idx][k] = 0 */
				mpq_set_ui(TN[idx][k], 0, 1);

				/* sum_r B[idx][0][r] * TN[idx_r][k] for n[r]>0 */
				for (r = 0; r < R; r++) {
					if (n[r] == 0) continue;
					n[r]--;
					int idx_r = popindex(n, R, planesizes);
					n[r]++;
					mpq_mul(tmp_q, B[idx][0][r], TN[idx_r][k]);
					mpq_add(TN[idx][k], TN[idx][k], tmp_q);
				}

				/* Add (1/alpha[k]) * B[idx][0][k] * (1 - Unp) */
				mpq_set_ui(tmp_q, 1, 1);
				mpq_sub(tmp_q, tmp_q, Unp); /* tmp_q = 1 - Unp */
				mpq_mul(tmp_q, tmp_q, B[idx][0][k]);
				mpq_set_z(tmp_q2, alpha[k]);
				mpq_div(tmp_q, tmp_q, tmp_q2); /* tmp_q = (1/alpha[k]) * B[idx][0][k] * (1-Unp) */
				mpq_add(TN[idx][k], TN[idx][k], tmp_q);

				mpq_clear(Unp);
			}
		}
	}
	while (!nextpop(n, qn->N, R));

	/* Clear progress line */
	fprintf(stderr, "\n");

	/* Extract final results at last population index */
	int final_idx = popindex(qn->N, R, planesizes);

	/* Q_out[m][r] = QN[final_idx][m][r] */
	for (s = 0; s < 2; s++) {
		for (r = 0; r < R; r++) {
			mpq_set(Q[s][r], QN[final_idx][s][r]);
		}
	}

	/* X[r] = TN[final_idx][r] */
	for (r = 0; r < R; r++) {
		mpq_set(X[r], TN[final_idx][r]);
	}

	/* B_out[m][r] = B[final_idx][m][r] (if B_out is not NULL) */
	if (B_out != NULL) {
		for (s = 0; s < 2; s++) {
			for (r = 0; r < R; r++) {
				mpq_set(B_out[s][r], B[final_idx][s][r]);
			}
		}
	}

	/* G = G_ca[final_idx] */
	mpq_set(G, G_ca[final_idx]);

	/* Compute logG */
	mpf_t G_mpf;
	mpf_init(G_mpf);
	mpf_set_q(G_mpf, G);
	double logG_final = log(mpf_get_d(G_mpf));
	mpf_clear(G_mpf);

	/* Free all memory */
	for (int i = 0; i < nstates; i++) {
		for (s = 0; s < 2; s++) {
			for (r = 0; r < R; r++) {
				mpq_clear(B[i][s][r]);
				mpq_clear(QN[i][s][r]);
			}
			free(B[i][s]);
			free(QN[i][s]);
		}
		free(B[i]);
		free(QN[i]);
		for (r = 0; r < R; r++)
			mpq_clear(TN[i][r]);
		free(TN[i]);
		mpq_clear(G_ca[i]);
		mpq_clear(V[i]);
	}
	free(B);
	free(QN);
	free(TN);
	free(G_ca);
	free(V);

	mpq_clear(tmp_q);
	mpq_clear(tmp_q2);
	mpq_clear(tmp_q3);
	mpq_clear(tmp_q4);
	mpz_clear(tmp_z);
	mpz_clear(tmp_z2);
	mpz_clear(scale_np);
	mpz_clear(scale_pr_k);
	mpq_clear(denom_q);

	free(n);
	free(planesizes);

	return logG_final;
}
