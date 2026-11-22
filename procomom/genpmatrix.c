#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "procomom.h"

/* Count the number of rows in the P matrices for class r */
static int count_rows(combsrep* Dn, int r, int R, int M)
{
	int W = Dn->stride;
	int numRows = 0;
	int d;

	for (d = 0; d < Dn->card; d++) {
		if (int_vecsubsum(Dn->combs[d], r-1, R-2) > 0) {
			/* Branch A: propagation, one row per component */
			numRows += W;
		} else {
			/* Branch B */
			if (int_vecsubsum(Dn->combs[d], 0, r-1) < M) {
				/* CE (one per station carried) + PC + extra PC */
				numRows += (W - 1) + (r - 1) + 1;
			} else {
				/* Boundary: extra PC only */
				numRows += 1;
			}
		}
	}

	return numRows;
}

/* genpmatrix - generate the 4 dense rectangular matrices A, B, DC, DD
   Port of MATLAB dmomprob.m nested function genpmatrix.
   A, B, DC, DD are [numRows x basisSize].
   Shift direction: UP (Dn[d] + e_s), matching CoMoM setupls.c.
   C and D matrices from MATLAB are always zero (nk=0), so omitted. */
PMatrices* genpmatrix(combsrep* Dn, qnmodel* qnm, int* Ncur, int r)
{
	int M = qnm->M;
	int R = qnm->R;
	mpz_t** L = qnm->L;
	mpz_t* Z = qnm->Z;

	int W = Dn->stride;            /* base + stations 1..W-1 */
	int basisSize = Dn->card * W;
	int numRows = count_rows(Dn, r, R, M);

	PMatrices* pm = (PMatrices*)calloc(1, sizeof(PMatrices));
	pm->numRows = numRows;
	pm->basisSize = basisSize;

	/* Allocate dense matrices (mpq_matzeros allocates rows+1 x cols+1) */
	pm->A = mpq_matzeros(numRows, basisSize);
	pm->B = mpq_matzeros(numRows, basisSize);
	pm->DC = mpq_matzeros(numRows, basisSize);
	pm->DD = mpq_matzeros(numRows, basisSize);

	int row = 0;
	int d, k, s;
	int* shifted = (int*)calloc(R, sizeof(int));
	mpz_t neg_val;
	mpz_init(neg_val);

	for (d = 0; d < Dn->card; d++) {
		if (int_vecsubsum(Dn->combs[d], r-1, R-2) > 0) {
			/*
			 * Branch A: propagation through class boundaries
			 * MATLAB: sum(Dn(d,r:R-1)) > 0
			 */
			for (k = 0; k < W; k++) {
				/* A[row, phash(Dn[d], k+1)] = 1 */
				int col_A = d * W + k;  /* = phash(Dn, Dn[d], k+1) */
				mpq_set_si(pm->A[row][col_A], 1, 1);

				if (int_vecsubsum(Dn->combs[d], r, R-2) > 0) {
					/* Classes beyond r have nonzero: B = A (identity copy) */
					mpq_set_si(pm->B[row][col_A], 1, 1);
				} else {
					/* Only class r has nonzero: B uses Dn[d] - e_r */
					int j;
					for (j = 0; j < R; j++) shifted[j] = Dn->combs[d][j];
					shifted[r-1]--;

					int col_B = phash(Dn, shifted, k+1);
					if (col_B >= 0 && col_B < basisSize)
						mpq_set_si(pm->B[row][col_B], 1, 1);
				}

				row++;
			}
		} else {
			/*
			 * Branch B: CE and PC equations
			 * sum(Dn(d,r:R-1)) == 0
			 */

			if (int_vecsubsum(Dn->combs[d], 0, r-1) < M) {
				/*
				 * CE equations: one per station carried, k=1..W-1
				 * component k against the base component
				 * Shift direction: UP (Dn[d] + e_s), matching CoMoM
				 */
				for (k = 1; k <= W-1; k++) {
					/* A: station k component = 1, base component = -1 */
					mpq_set_si(pm->A[row][d*W + k], 1, 1);
					mpq_set_si(pm->A[row][d*W + 0], -1, 1);

					/* Shifted terms for classes 1..r-1 */
					for (s = 1; s <= r-1; s++) {
						int j;
						for (j = 0; j < R; j++) shifted[j] = Dn->combs[d][j];
						shifted[s-1]++;  /* Dn[d] + e_s (UP shift, matching CoMoM) */

						int col = phash(Dn, shifted, k+1);
						if (col >= 0 && col < basisSize) {
							/* A: -L(k,s) at shifted, station k+1 */
							mpz_neg(neg_val, L[k-1][s-1]);
							mpq_set_z(pm->A[row][col], neg_val);
						}
					}

					/* B: L(k,r) at Dn[d], station k component */
					mpq_set_z(pm->B[row][d*W + k], L[k-1][r-1]);

					row++;
				}

				/*
				 * PC equations for classes s=1..r-1
				 * Population constraint with Z(s) and L(k,s)
				 * Shift direction: UP (Dn[d] + e_s), matching CoMoM
				 */
				for (s = 1; s <= r-1; s++) {
					/* nd(s) = Ncur[s-1] - Dn[d][s-1] */
					int nd_s = Ncur[s-1] - Dn->combs[d][s-1];

					/* A: nd(s) * base component of Dn[d] */
					mpq_set_si(pm->A[row][d*W + 0], nd_s, 1);

					/* Shifted = Dn[d] + e_s (UP shift, matching CoMoM) */
					int j;
					for (j = 0; j < R; j++) shifted[j] = Dn->combs[d][j];
					shifted[s-1]++;

					{
						int col_base = phash(Dn, shifted, 1);
						if (col_base >= 0 && col_base < basisSize) {
							/* A: -Z(s) at base component of shifted */
							mpz_neg(neg_val, Z[s-1]);
							mpq_set_z(pm->A[row][col_base], neg_val);

							/* DC: L(M,s) at base component of shifted */
							mpq_set_z(pm->DC[row][col_base], L[M-1][s-1]);
						}

						/* A: -m'(k)*L(k,s) for every station k=1..M at the
						 * station-k component of shifted.  The multiplicity is
						 * what comom/setupls.c:103 carries and this file was
						 * missing.  The station sum runs over the model with the
						 * reference COPY removed, so station M carries mi_M - 1
						 * (see the PC identity in promom/README.md); that term
						 * vanishes for a single-server reference, which is why
						 * it went unnoticed. */
						for (k = 1; k <= W-1; k++) {
							int col = phash(Dn, shifted, k+1);
							int mult = qnm->mi[k-1] - (k == M ? 1 : 0);
							if (col >= 0 && col < basisSize && mult != 0) {
								mpz_mul_si(neg_val, L[k-1][s-1], -mult);
								mpq_set_z(pm->A[row][col], neg_val);
							}
						}
					}

					/* B: zeroed (already 0 from mpq_matzeros) */

					row++;
				}
			}

			/*
			 * Extra PC for class r (always added in Branch B)
			 * This is the population constraint for the current class
			 */
			{
				/* nd(r) = Ncur[r-1] - Dn[d][r-1] */
				int nd_r = Ncur[r-1] - Dn->combs[d][r-1];

				/* A: nd(r) at base component of Dn[d] */
				mpq_set_si(pm->A[row][d*W + 0], nd_r, 1);

				/* B: Z(r) at base component of Dn[d] */
				mpq_set_z(pm->B[row][d*W + 0], Z[r-1]);

				/* B: m'(k)*L(k,r) for every station k=1..M at the station-k
				 * component of Dn[d], m'(M) = mi_M - 1 as above.
				 * Matches comom/setupls.c:135. */
				for (k = 1; k <= W-1; k++) {
					int mult = qnm->mi[k-1] - (k == M ? 1 : 0);
					if (mult == 0) continue;
					mpz_mul_si(neg_val, L[k-1][r-1], mult);
					mpq_set_z(pm->B[row][d*W + k], neg_val);
				}

				/* DD: L(M,r) at base component of Dn[d] */
				mpq_set_z(pm->DD[row][d*W + 0], L[M-1][r-1]);

				row++;
			}
		}
	}

	mpz_clear(neg_val);
	free(shifted);
	return pm;
}
