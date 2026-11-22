#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "mvaldmx.h"

/*
 * MVA-LD-MX solver: load-dependent MVA for mixed open/closed networks.
 * Port of pfqn_mvaldmx.m from LINE-dev.
 *
 * All arithmetic is exact rational (mpq_t).
 */

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) for k=0..Nt-1 */
void mvaldmx_auto_mu(qnmodel *qn)
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

/* hashpop: hash a population vector nvec[0..C-1] where 0<=nvec[c]<=Nc[c]
 * returns 0-based index = sum_c(prods[c]*nvec[c])
 * where prods[c] = prod_{j=0}^{c-1} (Nc[j]+1)
 * MATLAB returns 1-based, we return 0-based. */
static int hashpop(int *nvec, int *prods, int C)
{
	int idx = 0;
	for (int c = 0; c < C; c++)
		idx += prods[c] * nvec[c];
	return idx;
}

/* Advance population vector nvec (0<=nvec<=Nc) in pprod order.
 * Returns 0 if advanced, -1 if nvec == Nc (done). */
static int pprod_next(int *nvec, int *Nc, int C)
{
	/* Check if nvec == Nc */
	int all_equal = 1;
	for (int c = 0; c < C; c++) {
		if (nvec[c] != Nc[c]) {
			all_equal = 0;
			break;
		}
	}
	if (all_equal)
		return -1;

	int s = C - 1;
	while (s >= 0 && nvec[s] == Nc[s]) {
		nvec[s] = 0;
		s--;
	}
	if (s < 0)
		return -1;
	nvec[s]++;
	return 0;
}

double mvaldmx_solve(qnmodel *qn, mpq_t *X, mpq_t **Q)
{
	int M = qn->M;
	int R = qn->R;
	int Nt = qn->Nt;
	int ist, c, n, r;

	/* Identify closed and open classes */
	int numClosed = 0, numOpen = 0;
	int *closedIdx = getclosedclasses(qn, &numClosed);
	int *openIdx = getopenclasses(qn, &numOpen);
	int C = numClosed;

	if (C == 0) {
		fprintf(stderr, "mvaldmx: no closed classes found\n");
		if (openIdx) free(openIdx);
		return -1.0;
	}

	/* Build Nc, Dc, Zc for closed classes */
	int *Nc = (int *)calloc(C, sizeof(int));
	int Nct = 0;
	for (c = 0; c < C; c++) {
		Nc[c] = qn->N[closedIdx[c]];
		Nct += Nc[c];
	}

	mpz_t **Dc = (mpz_t **)malloc(M * sizeof(mpz_t *));
	for (ist = 0; ist < M; ist++) {
		Dc[ist] = (mpz_t *)malloc(C * sizeof(mpz_t));
		for (c = 0; c < C; c++) {
			mpz_init_set(Dc[ist][c], qn->L[ist][closedIdx[c]]);
		}
	}

	mpq_t *Zc = (mpq_t *)malloc(C * sizeof(mpq_t));
	for (c = 0; c < C; c++) {
		mpq_init(Zc[c]);
		mpq_set_z(Zc[c], qn->Z[closedIdx[c]]);
	}

	/* Ensure lambda exists (set to 0 for closed classes if not present) */
	mpq_t *lambda = (mpq_t *)malloc(R * sizeof(mpq_t));
	for (r = 0; r < R; r++) {
		mpq_init(lambda[r]);
		if (qn->hasOpen && qn->lambda != NULL)
			mpq_set(lambda[r], qn->lambda[r]);
		else
			mpq_set_ui(lambda[r], 0, 1);
	}

	/* Extend mu by one column: mu[:,Nt+1] = mu[:,Nt] (limited load dependence).
	 * Original mu has Nt columns (0..Nt-1). Extended has Nt+1 columns (0..Nt). */
	int mu_ext_cols = Nt + 1;
	mpq_t **mu_ext = (mpq_t **)malloc(M * sizeof(mpq_t *));
	for (ist = 0; ist < M; ist++) {
		mu_ext[ist] = (mpq_t *)malloc(mu_ext_cols * sizeof(mpq_t));
		for (int k = 0; k < Nt; k++) {
			mpq_init(mu_ext[ist][k]);
			mpq_set(mu_ext[ist][k], qn->mu[ist][k]);
		}
		/* Extra column: copy last value */
		mpq_init(mu_ext[ist][Nt]);
		mpq_set(mu_ext[ist][Nt], qn->mu[ist][Nt - 1]);
	}

	/* Compute EC, E, Eprime via effective capacity */
	/* EC[ist][n-1] for n=1..mu_ext_cols, E[ist][n] for n=0..mu_ext_cols, Eprime[ist][n] for n=0..mu_ext_cols */
	mpq_t **EC = (mpq_t **)malloc(M * sizeof(mpq_t *));
	mpq_t **E_arr = (mpq_t **)malloc(M * sizeof(mpq_t *));
	mpq_t **Eprime = (mpq_t **)malloc(M * sizeof(mpq_t *));
	for (ist = 0; ist < M; ist++) {
		EC[ist] = (mpq_t *)malloc(mu_ext_cols * sizeof(mpq_t));
		E_arr[ist] = (mpq_t *)malloc((mu_ext_cols + 1) * sizeof(mpq_t));
		Eprime[ist] = (mpq_t *)malloc((mu_ext_cols + 1) * sizeof(mpq_t));
		for (n = 0; n < mu_ext_cols; n++)
			mpq_init(EC[ist][n]);
		for (n = 0; n <= mu_ext_cols; n++) {
			mpq_init(E_arr[ist][n]);
			mpq_init(Eprime[ist][n]);
		}
	}

	mvaldmx_ec(EC, E_arr, Eprime, lambda, qn->L, mu_ext, M, R, Nt, mu_ext_cols);

	/* Compute prods for hashing: prods[c] = prod_{j=0}^{c-1} (Nc[j]+1) */
	int *prods = (int *)calloc(C, sizeof(int));
	prods[0] = 1;
	for (c = 1; c < C; c++)
		prods[c] = prods[c - 1] * (Nc[c - 1] + 1);

	int hashsize = 1;
	for (c = 0; c < C; c++)
		hashsize *= (Nc[c] + 1);

	/* Allocate Pc[M][Nct+1][hashsize] */
	/* Pc[ist][n][h] for n=0..Nct, h=0..hashsize-1 */
	mpq_t ***Pc = (mpq_t ***)malloc(M * sizeof(mpq_t **));
	for (ist = 0; ist < M; ist++) {
		Pc[ist] = (mpq_t **)malloc((Nct + 1) * sizeof(mpq_t *));
		for (n = 0; n <= Nct; n++) {
			Pc[ist][n] = (mpq_t *)malloc(hashsize * sizeof(mpq_t));
			for (int h = 0; h < hashsize; h++)
				mpq_init(Pc[ist][n][h]);
		}
	}

	/* Allocate w[M][C][hashsize] and x[C][hashsize] */
	mpq_t ***w = (mpq_t ***)malloc(M * sizeof(mpq_t **));
	for (ist = 0; ist < M; ist++) {
		w[ist] = (mpq_t **)malloc(C * sizeof(mpq_t *));
		for (c = 0; c < C; c++) {
			w[ist][c] = (mpq_t *)malloc(hashsize * sizeof(mpq_t));
			for (int h = 0; h < hashsize; h++)
				mpq_init(w[ist][c][h]);
		}
	}
	mpq_t **x = (mpq_t **)malloc(C * sizeof(mpq_t *));
	for (c = 0; c < C; c++) {
		x[c] = (mpq_t *)malloc(hashsize * sizeof(mpq_t));
		for (int h = 0; h < hashsize; h++)
			mpq_init(x[c][h]);
	}

	/* Initialize: nvec = (0,...,0), Pc[ist][0][hash(0,...,0)] = 1 for all ist */
	int *nvec = (int *)calloc(C, sizeof(int));
	int h0 = hashpop(nvec, prods, C);
	for (ist = 0; ist < M; ist++)
		mpq_set_ui(Pc[ist][0][h0], 1, 1);

	/* Temporary variables */
	mpq_t tmp, tmp2, sum_w;
	mpq_init(tmp);
	mpq_init(tmp2);
	mpq_init(sum_w);

	/* nvec_c: copy of nvec with one component decremented */
	int *nvec_c = (int *)malloc(C * sizeof(int));

	/* Population recursion */
	/* pprod iteration: start from (0,...,0) and advance to Nc */
	do {
		int hnvec = hashpop(nvec, prods, C);
		int nc = 0;
		for (c = 0; c < C; c++)
			nc += nvec[c];

		/* Compute w (mean residence times) */
		for (ist = 0; ist < M; ist++) {
			for (c = 0; c < C; c++) {
				if (nvec[c] > 0) {
					/* nvec_c = nvec with nvec[c]-- */
					memcpy(nvec_c, nvec, C * sizeof(int));
					nvec_c[c]--;
					int hnvec_c = hashpop(nvec_c, prods, C);

					for (n = 1; n <= nc; n++) {
						/* w[ist][c][hnvec] += Dc[ist][c] * n * EC[ist][n-1] * Pc[ist][n-1][hnvec_c] */
						mpq_set_z(tmp, Dc[ist][c]);
						mpq_set_ui(tmp2, n, 1);
						mpq_mul(tmp, tmp, tmp2);
						mpq_mul(tmp, tmp, EC[ist][n - 1]);
						mpq_mul(tmp, tmp, Pc[ist][n - 1][hnvec_c]);
						mpq_add(w[ist][c][hnvec], w[ist][c][hnvec], tmp);
					}
				}
			}
		}

		/* Compute throughput x */
		for (c = 0; c < C; c++) {
			/* x[c][hnvec] = nvec[c] / (Zc[c] + sum_ist(w[ist][c][hnvec])) */
			if (nvec[c] > 0) {
				mpq_set_ui(sum_w, 0, 1);
				for (ist = 0; ist < M; ist++)
					mpq_add(sum_w, sum_w, w[ist][c][hnvec]);
				mpq_add(sum_w, Zc[c], sum_w);
				mpq_set_ui(x[c][hnvec], nvec[c], 1);
				mpq_div(x[c][hnvec], x[c][hnvec], sum_w);
			}
			/* else x[c][hnvec] = 0 (already initialized) */
		}

		/* Update Pc */
		for (ist = 0; ist < M; ist++) {
			for (n = 1; n <= nc; n++) {
				for (c = 0; c < C; c++) {
					if (nvec[c] > 0) {
						memcpy(nvec_c, nvec, C * sizeof(int));
						nvec_c[c]--;
						int hnvec_c = hashpop(nvec_c, prods, C);
						/* Pc[ist][n][hnvec] += Dc[ist][c] * EC[ist][n-1] * x[c][hnvec] * Pc[ist][n-1][hnvec_c] */
						mpq_set_z(tmp, Dc[ist][c]);
						mpq_mul(tmp, tmp, EC[ist][n - 1]);
						mpq_mul(tmp, tmp, x[c][hnvec]);
						mpq_mul(tmp, tmp, Pc[ist][n - 1][hnvec_c]);
						mpq_add(Pc[ist][n][hnvec], Pc[ist][n][hnvec], tmp);
					}
				}
			}
			/* Pc[ist][0][hnvec] = max(0, 1 - sum(Pc[ist][1..nc][hnvec])) */
			mpq_set_ui(tmp, 0, 1);
			for (n = 1; n <= nc; n++)
				mpq_add(tmp, tmp, Pc[ist][n][hnvec]);
			mpq_set_ui(Pc[ist][0][hnvec], 1, 1);
			mpq_sub(Pc[ist][0][hnvec], Pc[ist][0][hnvec], tmp);
			/* If negative, set to 0 */
			if (mpq_sgn(Pc[ist][0][hnvec]) < 0)
				mpq_set_ui(Pc[ist][0][hnvec], 0, 1);
		}
	} while (pprod_next(nvec, Nc, C) == 0);

	/* Extract performance metrics at nvec = Nc */
	int hnvec_Nc = hashpop(Nc, prods, C);

	/* Throughputs for closed classes */
	for (c = 0; c < C; c++)
		mpq_set(X[closedIdx[c]], x[c][hnvec_Nc]);

	/* Utilization for closed classes: u[ist][c] */
	mpq_t **u = (mpq_t **)malloc(M * sizeof(mpq_t *));
	for (ist = 0; ist < M; ist++) {
		u[ist] = (mpq_t *)malloc(C * sizeof(mpq_t));
		for (c = 0; c < C; c++) {
			mpq_init(u[ist][c]);
			if (Nc[c] > 0) {
				memcpy(nvec_c, Nc, C * sizeof(int));
				nvec_c[c]--;
				int hnvec_c = hashpop(nvec_c, prods, C);
				for (n = 1; n <= Nct; n++) {
					/* u[ist][c] += Dc[ist][c] * x[c][hnvec_Nc] * Eprime[ist][n-1] / E_arr[ist][n-1] * Pc[ist][n-1][hnvec_c] */
					if (mpq_sgn(E_arr[ist][n - 1]) != 0) {
						mpq_set_z(tmp, Dc[ist][c]);
						mpq_mul(tmp, tmp, x[c][hnvec_Nc]);
						mpq_mul(tmp, tmp, Eprime[ist][n - 1]);
						mpq_div(tmp, tmp, E_arr[ist][n - 1]);
						mpq_mul(tmp, tmp, Pc[ist][n - 1][hnvec_c]);
						mpq_add(u[ist][c], u[ist][c], tmp);
					}
				}
			}
		}
	}

	/* Response times (CN) and queue lengths for closed classes */
	/* CN[ist][c] = w[ist][c][hnvec_Nc] */
	/* QN[ist][c] = X[closedIdx[c]] * CN[ist][c] */
	for (c = 0; c < C; c++) {
		for (ist = 0; ist < M; ist++) {
			mpq_mul(Q[ist][closedIdx[c]], X[closedIdx[c]], w[ist][c][hnvec_Nc]);
		}
	}

	/* Open class metrics */
	for (int oi = 0; oi < numOpen; oi++) {
		r = openIdx[oi];
		mpq_set(X[r], lambda[r]);
		for (ist = 0; ist < M; ist++) {
			/* QN[ist][r] = sum_{n=0}^{Nct} lambda[r] * D[ist][r] * (n+1) * EC[ist][n+1-1] * Pc[ist][n][hnvec_Nc] */
			/* EC[ist][n] for n=1..Nt_ec, so EC[ist][n+1-1] = EC[ist][n] */
			mpq_set_ui(Q[ist][r], 0, 1);
			for (n = 0; n <= Nct; n++) {
				/* Need EC[ist][n+1-1] = EC[ist][n], which requires n+1 <= mu_ext_cols */
				if (n + 1 <= mu_ext_cols) {
					mpq_set(tmp, lambda[r]);
					mpq_set_z(tmp2, qn->L[ist][r]);
					mpq_mul(tmp, tmp, tmp2);
					mpq_set_ui(tmp2, n + 1, 1);
					mpq_mul(tmp, tmp, tmp2);
					mpq_mul(tmp, tmp, EC[ist][n]); /* EC[ist][n] = EC for rate level n+1 */
					mpq_mul(tmp, tmp, Pc[ist][n][hnvec_Nc]);
					mpq_add(Q[ist][r], Q[ist][r], tmp);
				}
			}
		}
	}

	/* Compute logG: not directly computed by MVA, return 0.0 */
	double logG = 0.0;

	/* Cleanup */
	for (ist = 0; ist < M; ist++) {
		for (c = 0; c < C; c++) {
			mpq_clear(u[ist][c]);
		}
		free(u[ist]);
	}
	free(u);

	for (ist = 0; ist < M; ist++) {
		for (n = 0; n <= Nct; n++) {
			for (int h = 0; h < hashsize; h++)
				mpq_clear(Pc[ist][n][h]);
			free(Pc[ist][n]);
		}
		free(Pc[ist]);
	}
	free(Pc);

	for (ist = 0; ist < M; ist++) {
		for (c = 0; c < C; c++) {
			for (int h = 0; h < hashsize; h++)
				mpq_clear(w[ist][c][h]);
			free(w[ist][c]);
		}
		free(w[ist]);
	}
	free(w);

	for (c = 0; c < C; c++) {
		for (int h = 0; h < hashsize; h++)
			mpq_clear(x[c][h]);
		free(x[c]);
	}
	free(x);

	for (ist = 0; ist < M; ist++) {
		for (n = 0; n < mu_ext_cols; n++)
			mpq_clear(EC[ist][n]);
		for (n = 0; n <= mu_ext_cols; n++) {
			mpq_clear(E_arr[ist][n]);
			mpq_clear(Eprime[ist][n]);
		}
		free(EC[ist]);
		free(E_arr[ist]);
		free(Eprime[ist]);
	}
	free(EC);
	free(E_arr);
	free(Eprime);

	for (ist = 0; ist < M; ist++) {
		for (int k = 0; k < mu_ext_cols; k++)
			mpq_clear(mu_ext[ist][k]);
		free(mu_ext[ist]);
	}
	free(mu_ext);

	for (ist = 0; ist < M; ist++) {
		for (c = 0; c < C; c++)
			mpz_clear(Dc[ist][c]);
		free(Dc[ist]);
	}
	free(Dc);

	for (c = 0; c < C; c++)
		mpq_clear(Zc[c]);
	free(Zc);

	for (r = 0; r < R; r++)
		mpq_clear(lambda[r]);
	free(lambda);

	free(Nc);
	free(prods);
	free(nvec);
	free(nvec_c);
	if (closedIdx) free(closedIdx);
	if (openIdx) free(openIdx);

	mpq_clear(tmp);
	mpq_clear(tmp2);
	mpq_clear(sum_w);

	return logG;
}
