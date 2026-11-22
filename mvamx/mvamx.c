#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "mvamx.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

/**
 * Closed MVA with mpq_t demands.
 *
 * Runs the standard multiclass MVA population recursion on a closed sub-model
 * whose service demands Dc[m][c] are mpq_t rationals (already scaled by
 * 1/(1 - UNt[m])).
 *
 * @param M       Number of stations
 * @param Rc      Number of closed classes
 * @param Nc      Population vector for closed classes [Rc]
 * @param Zc      Think times for closed classes [Rc] (mpz_t)
 * @param Dc      Scaled demands [M][Rc] (mpq_t)
 * @param mi      Station multiplicities [M]
 * @param Xc      Output: throughputs for closed classes [Rc] (mpq_t, pre-initialized)
 * @param Qc      Output: queue lengths [M][Rc] (mpq_t, pre-initialized)
 * @param Gc      Output: normalizing constant (mpq_t, pre-initialized)
 * @return        log(G) of the closed sub-model
 */
static double closed_mva_mpq(int M, int Rc, int *Nc, mpz_t *Zc,
                              mpq_t **Dc, int *mi,
                              mpq_t *Xc, mpq_t **Qc, mpq_t Gc)
{
	struct rusage ruse;
	double t_start = CPUTIME;
	int r, m, t;

	int *n = (int*)initpop(Rc);
	if (n == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory in closed_mva_mpq\n");
		return -1.0;
	}
	int *planesizes = getplanesizes(Nc, Rc);
	if (planesizes == NULL) {
		free(n);
		fprintf(stderr, "Error: Failed to allocate memory in closed_mva_mpq\n");
		return -1.0;
	}

	int planesize0 = planesizes[0];

	/* q[popindex][m] = total queue length at station m for population n */
	mpq_t **q = (mpq_t**)calloc(planesize0, sizeof(mpq_t*));
	if (q == NULL) {
		free(n);
		free(planesizes);
		fprintf(stderr, "Error: Failed to allocate memory for q array (%d elements). Model may be too large.\n", planesize0);
		return -1.0;
	}
	for (t = 0; t < planesize0; t++) {
		q[t] = (mpq_t*)calloc(M, sizeof(mpq_t));
		if (q[t] == NULL) {
			for (int i = 0; i < t; i++) {
				for (m = 0; m < M; m++)
					mpq_clear(q[i][m]);
				free(q[i]);
			}
			free(q);
			free(n);
			free(planesizes);
			fprintf(stderr, "Error: Failed to allocate memory for q[%d]. Model may be too large.\n", t);
			return -1.0;
		}
		for (m = 0; m < M; m++)
			mpq_init(q[t][m]);
	}

	int curindex = -1;

	/* G accumulator: G starts at 1 at population (0,...,0) */
	mpq_t G_accumulator;
	mpq_init(G_accumulator);
	mpq_set_ui(G_accumulator, 1, 1);

	do {
		curindex++;

		/* Progress output */
		int *widths = (int*)calloc(Rc, sizeof(int));
		for (r = 0; r < Rc; r++) {
			int val = Nc[r];
			widths[r] = 1;
			while (val >= 10) { widths[r]++; val /= 10; }
		}
		fprintf(stderr, "\rn=(");
		for (r = 0; r < Rc; r++) {
			fprintf(stderr, "%*d", widths[r], n[r]);
			if (r < Rc - 1) fprintf(stderr, ",");
		}
		double t_current = CPUTIME;
		fprintf(stderr, ") - Time: %.2f s  ", t_current - t_start);
		fflush(stderr);
		free(widths);

		for (r = 0; r < Rc; r++) {
			if (n[r] != 0) {
				n[r]--;
				int index_1r = popindex(n, Rc, planesizes);
				n[r]++;

				/* Cr = Zc[r] + sum_m Dc[m][r] * (mi[m] + q[index_1r][m]) */
				mpq_t Cr;
				mpq_init(Cr);
				mpq_set_z(Cr, Zc[r]); /* Cr = Z[r] (integer -> rational) */

				for (m = 0; m < M; m++) {
					/* mi_plus_q = mi[m] + q[index_1r][m] */
					mpq_t mi_plus_q;
					mpq_init(mi_plus_q);
					mpq_set_ui(mi_plus_q, mi[m], 1);
					mpq_add(mi_plus_q, mi_plus_q, q[index_1r][m]);

					/* term = Dc[m][r] * mi_plus_q */
					mpq_t term;
					mpq_init(term);
					mpq_mul(term, Dc[m][r], mi_plus_q);
					mpq_add(Cr, Cr, term);

					mpq_clear(term);
					mpq_clear(mi_plus_q);
				}

				/* Xc[r] = n[r] / Cr */
				mpq_set_ui(Xc[r], n[r], 1);
				mpq_div(Xc[r], Xc[r], Cr);

				/* Update queue lengths */
				for (m = 0; m < M; m++) {
					mpq_t mi_plus_q;
					mpq_init(mi_plus_q);
					mpq_set_ui(mi_plus_q, mi[m], 1);
					mpq_add(mi_plus_q, mi_plus_q, q[index_1r][m]);

					mpq_t update;
					mpq_init(update);
					mpq_mul(update, Dc[m][r], mi_plus_q);
					mpq_mul(update, update, Xc[r]);
					mpq_add(q[curindex][m], q[curindex][m], update);

					mpq_clear(update);
					mpq_clear(mi_plus_q);
				}

				mpq_clear(Cr);
			}
		}

		/* Accumulate G using the same logic as mva-multi.c */
		int nonzero_count = 0;
		int last_nnz = -1;
		for (r = 0; r < Rc; r++) {
			if (n[r] > 0) {
				last_nnz = r;
				nonzero_count++;
			}
		}
		if (nonzero_count > 0 && last_nnz >= 0) {
			int sumn = 0, sumN = 0;
			for (r = 0; r < last_nnz; r++) {
				sumn += n[r];
				sumN += Nc[r];
			}
			int sumnprime = 0;
			for (r = last_nnz + 1; r < Rc; r++)
				sumnprime += n[r];

			if (sumn == sumN && sumnprime == 0) {
				mpq_div(G_accumulator, G_accumulator, Xc[last_nnz]);
			}
		}

	} while (!nextpop(n, Nc, Rc));

	fprintf(stderr, "\n");

	mpq_set(Gc, G_accumulator);
	mpq_clear(G_accumulator);

	/* Compute final queue lengths at the full population Nc */
	for (m = 0; m < M; m++) {
		for (r = 0; r < Rc; r++) {
			int *n_minus = (int*)malloc(Rc * sizeof(int));
			for (int i = 0; i < Rc; i++)
				n_minus[i] = Nc[i];
			n_minus[r]--;

			int idx = popindex(n_minus, Rc, planesizes);

			mpq_t mi_plus_q;
			mpq_init(mi_plus_q);
			mpq_set_ui(mi_plus_q, mi[m], 1);
			mpq_add(mi_plus_q, mi_plus_q, q[idx][m]);

			mpq_mul(Qc[m][r], Dc[m][r], mi_plus_q);
			mpq_mul(Qc[m][r], Qc[m][r], Xc[r]);

			mpq_clear(mi_plus_q);
			free(n_minus);
		}
	}

	/* Free memory */
	for (t = 0; t < planesize0; t++) {
		for (m = 0; m < M; m++)
			mpq_clear(q[t][m]);
		free(q[t]);
	}
	free(q);
	free(n);
	free(planesizes);

	/* Compute logG from final G */
	mpf_t G_mpf;
	mpf_init(G_mpf);
	mpf_set_q(G_mpf, Gc);
	double logG = log(mpf_get_d(G_mpf));
	mpf_clear(G_mpf);

	return logG;
}

double mvamx(qnmodel *qn, mpq_t *X, mpq_t **Q, mpq_t G)
{
	int m, r;
	int M = qn->M;

	/* ---- Identify open and closed classes ---- */
	int numOpen = 0, numClosed = 0;
	int *openIdx = getopenclasses(qn, &numOpen);
	int *closedIdx = getclosedclasses(qn, &numClosed);

	/* ---- Step 1: Open-class metrics ---- */
	/* X[r] = lambda[r],  U[m][r] = lambda[r] * L[m][r] */
	/* UNt[m] = sum over open classes of U[m][r]           */
	mpq_t *UNt = (mpq_t*)malloc(M * sizeof(mpq_t));
	for (m = 0; m < M; m++) {
		mpq_init(UNt[m]);
		mpq_set_ui(UNt[m], 0, 1);
	}

	for (int i = 0; i < numOpen; i++) {
		r = openIdx[i];
		mpq_set(X[r], qn->lambda[r]); /* X[r] = lambda[r] */
		for (m = 0; m < M; m++) {
			/* U = lambda[r] * L[m][r] */
			mpq_t U;
			mpq_init(U);
			mpq_set_z(U, qn->L[m][r]);
			mpq_mul(U, U, qn->lambda[r]);
			mpq_add(UNt[m], UNt[m], U);
			mpq_clear(U);
		}
	}

	/* Check stability: UNt[m] < 1 for all stations */
	for (m = 0; m < M; m++) {
		if (mpq_cmp_ui(UNt[m], 1, 1) >= 0) {
			fprintf(stderr, "Error: Station %d is unstable (utilization >= 1)\n", m + 1);
			for (m = 0; m < M; m++) mpq_clear(UNt[m]);
			free(UNt);
			if (openIdx) free(openIdx);
			if (closedIdx) free(closedIdx);
			return -1.0;
		}
	}

	/* ---- Step 2: Solve closed sub-model ---- */
	double logG = 0.0;
	mpq_set_ui(G, 1, 1); /* default G=1 if no closed classes */

	/* We will store Qc[m][c] for the closed sub-model so we can compute
	   open-class cycle times afterwards. */
	mpq_t **Qc = NULL; /* [M][numClosed] */

	if (numClosed > 0) {
		/* Build Nc, Zc, Dc for the closed sub-model */
		int *Nc = (int*)malloc(numClosed * sizeof(int));
		mpz_t *Zc = (mpz_t*)malloc(numClosed * sizeof(mpz_t));
		mpq_t **Dc = (mpq_t**)malloc(M * sizeof(mpq_t*));

		for (int c = 0; c < numClosed; c++) {
			Nc[c] = qn->N[closedIdx[c]];
			mpz_init_set(Zc[c], qn->Z[closedIdx[c]]);
		}

		/* Dc[m][c] = L[m][closedIdx[c]] / (1 - UNt[m]) */
		for (m = 0; m < M; m++) {
			Dc[m] = (mpq_t*)malloc(numClosed * sizeof(mpq_t));
			for (int c = 0; c < numClosed; c++) {
				mpq_init(Dc[m][c]);
				/* numerator = L[m][closedIdx[c]] */
				mpq_set_z(Dc[m][c], qn->L[m][closedIdx[c]]);
				/* denominator factor = (1 - UNt[m]) */
				mpq_t one_minus_u;
				mpq_init(one_minus_u);
				mpq_set_ui(one_minus_u, 1, 1);
				mpq_sub(one_minus_u, one_minus_u, UNt[m]);
				/* Dc[m][c] = L / (1 - UNt) */
				mpq_div(Dc[m][c], Dc[m][c], one_minus_u);
				mpq_clear(one_minus_u);
			}
		}

		/* Allocate Xc, Qc */
		mpq_t *Xc = (mpq_t*)malloc(numClosed * sizeof(mpq_t));
		Qc = (mpq_t**)malloc(M * sizeof(mpq_t*));
		for (int c = 0; c < numClosed; c++)
			mpq_init(Xc[c]);
		for (m = 0; m < M; m++) {
			Qc[m] = (mpq_t*)malloc(numClosed * sizeof(mpq_t));
			for (int c = 0; c < numClosed; c++)
				mpq_init(Qc[m][c]);
		}

		/* Run closed MVA with mpq_t demands */
		logG = closed_mva_mpq(M, numClosed, Nc, Zc, Dc, qn->mi, Xc, Qc, G);

		/* Map closed-class results back to full arrays */
		for (int c = 0; c < numClosed; c++) {
			mpq_set(X[closedIdx[c]], Xc[c]);
			for (m = 0; m < M; m++)
				mpq_set(Q[m][closedIdx[c]], Qc[m][c]);
		}

		/* Free closed sub-model temporaries */
		for (int c = 0; c < numClosed; c++) {
			mpq_clear(Xc[c]);
			mpz_clear(Zc[c]);
		}
		for (m = 0; m < M; m++) {
			for (int c = 0; c < numClosed; c++)
				mpq_clear(Dc[m][c]);
			free(Dc[m]);
		}
		free(Xc);
		free(Dc);
		free(Nc);
		free(Zc);
	}

	/* ---- Step 3: Open-class cycle times and queue lengths ---- */
	/* CN[m][r] = L[m][r] * (1 + sum_c Qc[m][c]) / (1 - UNt[m])
	   Q[m][r]  = CN[m][r] * lambda[r]                              */
	for (int i = 0; i < numOpen; i++) {
		r = openIdx[i];
		for (m = 0; m < M; m++) {
			/* one_plus_sumQc = 1 + sum_c Qc[m][c] */
			mpq_t one_plus_sumQc;
			mpq_init(one_plus_sumQc);
			mpq_set_ui(one_plus_sumQc, 1, 1);

			if (Qc != NULL) {
				for (int c = 0; c < numClosed; c++) {
					mpq_add(one_plus_sumQc, one_plus_sumQc, Qc[m][c]);
				}
			}

			/* one_minus_u = 1 - UNt[m] */
			mpq_t one_minus_u;
			mpq_init(one_minus_u);
			mpq_set_ui(one_minus_u, 1, 1);
			mpq_sub(one_minus_u, one_minus_u, UNt[m]);

			/* CN = L[m][r] * one_plus_sumQc / one_minus_u */
			mpq_t CN;
			mpq_init(CN);
			mpq_set_z(CN, qn->L[m][r]);
			mpq_mul(CN, CN, one_plus_sumQc);
			mpq_div(CN, CN, one_minus_u);

			/* Q[m][r] = CN * lambda[r] */
			mpq_mul(Q[m][r], CN, qn->lambda[r]);

			mpq_clear(CN);
			mpq_clear(one_minus_u);
			mpq_clear(one_plus_sumQc);
		}
	}

	/* ---- Cleanup ---- */
	if (Qc != NULL) {
		/* Qc was already cleared for closed-class entries above in Dc cleanup?
		   No -- Qc is separate from Dc. We need to free Qc here. */
		for (m = 0; m < M; m++) {
			/* Qc[m][c] entries are still live -- but we already mapped them to
			   Q[][] above. We need to clear and free. BUT: the entries were
			   cleared by the closed-class block only for Dc and Xc, not Qc.
			   We must clear them here. */
			for (int c = 0; c < numClosed; c++)
				mpq_clear(Qc[m][c]);
			free(Qc[m]);
		}
		free(Qc);
	}

	for (m = 0; m < M; m++)
		mpq_clear(UNt[m]);
	free(UNt);
	if (openIdx) free(openIdx);
	if (closedIdx) free(closedIdx);

	return logG;
}
