#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

/* Port of setup1.m: build the six matrices A,B,C,D,E,F for the prefix
 * sub-model on queues 1..m (queue m removed to reach level m-1), for the
 * first r classes.  L is m x r (mpz demands), N and Z length r.
 *
 * Indexing note: MATLAB is 1-based and uses (i-1)*R+1+s with s=0..r-1;
 * the C form is i*r + s with 0-based combination index i.  The row count
 * of the CE+PC system is overdetermined and is computed in a first pass.
 */

static void set_q_from_z(mpq_t dst, mpz_t z)      { mpq_set_z(dst, z); }
static void set_q_neg_z(mpq_t dst, mpz_t z)       { mpq_set_z(dst, z); mpq_neg(dst, dst); }

LS1* setup1(mpz_t** L, int m, int r, int* N, mpz_t* Z)
{
	int i, j, s, row, t;
	LS1* ls = (LS1*) calloc(1, sizeof(LS1));

	int** Ik = sortbynnzpos(multichoose(m, r-1),   nck(m+(r-1)-1, m-1),   m);
	int** Ii = sortbynnzpos(multichoose(m-1, r-1), nck((m-1)+(r-1)-1, m-2), m-1);
	int** I  = sortbynnzpos(multichoose(m, r-2),   nck(m+(r-2)-1, m-1),   m);
	int szIk = nck(m+(r-1)-1, m-1);
	int szIi = nck((m-1)+(r-1)-1, m-2);
	int szI  = nck(m+(r-2)-1, m-1);
	int colA = szIk * r, colC = szIi * r, colD = szI * r;

	ls->szIk = szIk; ls->szIi = szIi; ls->szI = szI;
	ls->colA = colA; ls->colC = colC; ls->colD = colD;

	/* Pass 1: count overdetermined rows of A,B,C,D. */
	int nrows = 0;
	for (i = 0; i < szIk; i++) if (Ik[i][m-1] == 0) nrows++;       /* subtractive, j=M */
	for (i = 0; i < szIk; i++) for (j = 0; j < m; j++) if (Ik[i][j] > 0) nrows++; /* additive */
	nrows += szI * (r - 1);                                        /* PC classes 1..r-1 */
	ls->nrows = nrows;

	ls->A = mpq_matzeros(nrows, colA);
	ls->B = mpq_matzeros(nrows, colA);
	ls->C = mpq_matzeros(nrows, colC);
	ls->D = mpq_matzeros(nrows, colD);
	ls->E = mpq_matzeros(colD, colA);   /* szI*r rows x szIk*r cols */
	ls->F = mpq_matzeros(colD, colD);

	row = 0;

	/* SUBTRACTIVE CONVOLUTION (queue m removed) */
	for (i = 0; i < szIk; i++) {
		if (Ik[i][m-1] == 0) {
			mpq_set_si(ls->A[row][i*r + 0], 1, 1);
			set_q_from_z(ls->B[row][i*r + 0], L[m-1][r-1]);          /* L(M,r) */
			t = (int) int_matmatchrow(Ii, szIi, m-1, Ik[i]);         /* Ik(i,1:M-1) */
			mpq_set_si(ls->C[row][t*r + 0], 1, 1);
			for (s = 1; s < r; s++) set_q_neg_z(ls->A[row][i*r + s], L[m-1][s-1]);
			row++;
		}
	}

	/* ADDITIVE CONVOLUTION */
	for (i = 0; i < szIk; i++) {
		for (j = 0; j < m; j++) {
			if (Ik[i][j] > 0) {
				mpq_set_si(ls->A[row][i*r + 0], 1, 1);
				set_q_from_z(ls->B[row][i*r + 0], L[j][r-1]);
				Ik[i][j]--;
				t = (int) int_matmatchrow(I, szI, m, Ik[i]);
				Ik[i][j]++;
				mpq_set_si(ls->D[row][t*r + 0], 1, 1);
				for (s = 1; s < r; s++) set_q_neg_z(ls->A[row][i*r + s], L[j][s-1]);
				row++;
			}
		}
	}

	/* POPULATION CONSTRAINTS, classes 1..r-1 */
	for (i = 0; i < szI; i++) {
		for (s = 1; s < r; s++) {
			mpq_set_si(ls->D[row][i*r + 0], -N[s-1], 1);            /* -N(s) */
			set_q_from_z(ls->D[row][i*r + s], Z[s-1]);              /* Z(s) */
			for (j = 0; j < m; j++) {
				I[i][j]++;
				t = (int) int_matmatchrow(Ik, szIk, m, I[i]);
				I[i][j]--;
				/* -(1 + I(i,j)) * L(j,s) */
				mpz_t v; mpz_init(v);
				mpz_mul_ui(v, L[j][s-1], (unsigned long)(1 + I[i][j]));
				mpz_neg(v, v);
				mpq_set_z(ls->A[row][t*r + s], v);
				mpz_clear(v);
			}
			row++;
		}
	}

	/* POPULATION CONSTRAINTS, class R (=r): E, F */
	row = 0;
	for (i = 0; i < szI; i++) {
		for (s = 0; s < r; s++) {
			set_q_from_z(ls->F[row][i*r + s], Z[r-1]);              /* Z(r) */
			for (j = 0; j < m; j++) {
				I[i][j]++;
				t = (int) int_matmatchrow(Ik, szIk, m, I[i]);
				I[i][j]--;
				mpz_t v; mpz_init(v);
				mpz_mul_ui(v, L[j][r-1], (unsigned long)(1 + I[i][j]));
				mpq_set_z(ls->E[row][t*r + s], v);
				mpz_clear(v);
			}
			row++;
		}
	}

	free(Ik); free(Ii); free(I);
	return ls;
}

void ls1_free(LS1* ls, int r)
{
	if (!ls) return;
	free(ls->A); free(ls->B); free(ls->C); free(ls->D); free(ls->E); free(ls->F);
	free(ls);
}
