#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "promom.h"

/* Combination helpers ---------------------------------------------------- */

qcombs* qcombs_new(int nvec, int level)
{
	qcombs* c = (qcombs*) calloc(1, sizeof(qcombs));
	c->n = nvec;
	c->k = level;
	c->card = nck(nvec + level - 1, level);
	c->combs = (int**) sortbynnzpos((int**) multichoose(nvec, level), c->card, nvec);
	return c;
}

void qcombs_free(qcombs* c)
{
	if (!c) return;
	free(c->combs);
	free(c);
}

int qcombs_index(qcombs* c, int* comb)
{
	long int pos = int_matmatchrow(c->combs, c->card, c->n, comb);
	return (int) pos;
}

/* Row count for the step system ------------------------------------------ */

static int count_rows(qcombs* Ik, qcombs* I, int M, int r)
{
	int w, j, ce = 0;
	for (w = 0; w < Ik->card; w++)
		for (j = 0; j < M; j++)
			if (Ik->combs[w][j] > 0) ce++;
	/* CE, one per (replica combination, occupied station);
	 * PC for classes s = 1..r-1, one per level-(r-1) combination;
	 * PC for class r, one per level-(r-1) combination and decrement. */
	return ce + I->card * (r - 1) + I->card * r;
}

/* genqmatrix -------------------------------------------------------------
 *
 * Assembles A, B, DC, DD for one population step of class r at population
 * Ncur, over the MoM basis described in promom.h.  The reference station is
 * station M (the caller rotates the model so that this is always true), so
 * the n-coupled entries carry L[M-1][s].
 *
 * With the n-coupled blocks dropped this is mom's own step system, laid out
 * densely instead of in BTF blocks: rows 1 are mom's CE rows
 * (mom/setupls.c:83-107), rows 2 its class-s population constraints
 * (setupls.c:192-202), rows 3 its class-r population constraints
 * (setupls.c:120-191, the B2r block).
 */
QMatrices* genqmatrix(qcombs* Ik, qcombs* I, qnmodel* qnm, int* Ncur, int r)
{
	int M = qnm->M;
	mpz_t** L = qnm->L;
	mpz_t*  Z = qnm->Z;
	int*    mi = qnm->mi;

	int cardk = Ik->card, cardi = I->card;
	int off = cardk * r;                 /* start of the level-(r-1) block */
	int ncols = (cardk + cardi) * r;
	int nrows = count_rows(Ik, I, M, r);

	QMatrices* qm = (QMatrices*) calloc(1, sizeof(QMatrices));
	qm->nrows = nrows;
	qm->ncols = ncols;
	qm->A  = mpq_matzeros(nrows, ncols);
	qm->B  = mpq_matzeros(nrows, ncols);
	qm->DC = mpq_matzeros(nrows, ncols);
	qm->DD = mpq_matzeros(nrows, ncols);

	int row = 0, w, i, j, s;
	mpz_t neg;
	mpz_init(neg);
	int* shifted = (int*) calloc(M, sizeof(int));

	/* ---- rows 1: convolution equations -------------------------------
	 * q(w, P) - sum_{s=1}^{r-1} L_js q(w, P-1_s) - q(w-1_j, P)
	 *     = L_jr q(w, P-1_r)
	 * The right-hand side is the previous point's level-r base component.
	 * No n-coupling: a replica of station j leaves the reference count
	 * untouched.                                                          */
	for (w = 0; w < cardk; w++) {
		for (j = 0; j < M; j++) {
			if (Ik->combs[w][j] == 0) continue;
			mpq_set_si(qm->A[row][w*r + 0], 1, 1);
			for (s = 1; s <= r-1; s++) {
				mpz_neg(neg, L[j][s-1]);
				mpq_set_z(qm->A[row][w*r + s], neg);
			}
			for (i = 0; i < M; i++) shifted[i] = Ik->combs[w][i];
			shifted[j]--;
			{
				int iw = qcombs_index(I, shifted);
				if (iw >= 0) mpq_set_si(qm->A[row][off + iw*r + 0], -1, 1);
			}
			mpq_set_z(qm->B[row][w*r + 0], L[j][r-1]);
			row++;
		}
	}

	/* ---- rows 2: population constraints for classes s = 1..r-1 -------
	 * P_s q(i, P) - Z_s q(i, P-1_s)
	 *     - sum_j (mi_j + i_j - [j=M]) L_js q(i+1_j, P-1_s)
	 *     = n L_Ms q(i, P-1_s, n-1)
	 * every term at the current point, so the n-coupling sits in DC.     */
	for (i = 0; i < cardi; i++) {
		for (s = 1; s <= r-1; s++) {
			mpq_set_si(qm->A[row][off + i*r + 0], Ncur[s-1], 1);
			mpz_neg(neg, Z[s-1]);
			mpq_set_z(qm->A[row][off + i*r + s], neg);
			for (j = 0; j < M; j++) {
				int t;
				for (t = 0; t < M; t++) shifted[t] = I->combs[i][t];
				shifted[j]++;
				{
					int iw = qcombs_index(Ik, shifted);
					int mult = mi[j] + I->combs[i][j] - (j == M-1 ? 1 : 0);
					if (iw < 0 || mult == 0) continue;
					mpz_mul_ui(neg, L[j][s-1], (unsigned long) mult);
					mpz_neg(neg, neg);
					mpq_set_z(qm->A[row][iw*r + s], neg);
				}
			}
			mpq_set_z(qm->DC[row][off + i*r + s], L[M-1][s-1]);
			row++;
		}
	}

	/* ---- rows 3: population constraint for class r -------------------
	 * applied at P - 1_s for each decrement s = 0..r-1:
	 * P_r q(i, P-1_s) = Z_r q(i, P-1_s-1_r)
	 *     + sum_j (mi_j + i_j - [j=M]) L_jr q(i+1_j, P-1_s-1_r)
	 *     + n L_Mr q(i, P-1_s-1_r, n-1)
	 * P - 1_s - 1_r is the previous point decremented by s, so the
	 * right-hand side reads the previous vector: B and DD.               */
	for (i = 0; i < cardi; i++) {
		for (s = 0; s <= r-1; s++) {
			mpq_set_si(qm->A[row][off + i*r + s], Ncur[r-1], 1);
			mpq_set_z(qm->B[row][off + i*r + s], Z[r-1]);
			for (j = 0; j < M; j++) {
				int t;
				for (t = 0; t < M; t++) shifted[t] = I->combs[i][t];
				shifted[j]++;
				{
					int iw = qcombs_index(Ik, shifted);
					int mult = mi[j] + I->combs[i][j] - (j == M-1 ? 1 : 0);
					if (iw < 0 || mult == 0) continue;
					mpz_mul_ui(neg, L[j][r-1], (unsigned long) mult);
					mpq_set_z(qm->B[row][iw*r + s], neg);
				}
			}
			mpq_set_z(qm->DD[row][off + i*r + s], L[M-1][r-1]);
			row++;
		}
	}

	if (row != nrows) {
		fprintf(stderr, "promom: internal error, filled %d of %d rows\n", row, nrows);
		exit(2);
	}

	mpz_clear(neg);
	free(shifted);
	return qm;
}

void free_qmatrices(QMatrices* qm)
{
	int i, j;
	for (i = 0; i <= qm->nrows; i++) {
		for (j = 0; j <= qm->ncols; j++) {
			mpq_clear(qm->A[i][j]);  mpq_clear(qm->B[i][j]);
			mpq_clear(qm->DC[i][j]); mpq_clear(qm->DD[i][j]);
		}
		free(qm->A[i]); free(qm->B[i]); free(qm->DC[i]); free(qm->DD[i]);
	}
	free(qm->A); free(qm->B); free(qm->DC); free(qm->DD);
	free(qm);
}
