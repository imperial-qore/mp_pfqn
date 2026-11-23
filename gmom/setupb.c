#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

/* Level system of the generalized MoM with variable branching factor
 * (Casale, Perf. Eval. 2011, Sec. 4.1.2 and 4.2), for one sub-model.
 *
 * The sub-model holds the d original queues q[0..d-1], each with unit
 * multiplicity, and carries the bases of levels l (unknowns, vl) and l-1
 * (vl1).  An entry of the level-l basis is the r-vector
 *   G(m'+delta, n), G(m'+delta, n-1_1), ..., G(m'+delta, n-1_{r-1}),
 * delta a multiset of l replicas over the d queues, combinations in the
 * sortbynnzpos(multichoose(d, l)) order.  The rows are, in order:
 *
 *   GCE  for each tagged queue (the last ntag of q) and each delta with
 *        delta_j = 0, the general convolution expression (20) removing
 *        queue j:
 *          G(m'+delta, n) - sum_{s<r} L_js G(m'+delta, n-1_s)
 *             = G(m'-1_j+delta, n) + L_jr G(m'+delta, n-1_r);
 *        the first right-hand term is component 0 of the child model's
 *        level-l basis, in its vl slot when the child carries level l
 *        and in its vl1 slot when it carries level l+1 (Theorem 2);
 *   CE   the same with a replica of queue j (delta_j > 0) removed, whose
 *        first right-hand term is in this model's level l-1 (D);
 *   PC   classes 1..r-1, as in setup1.
 *
 * E and F give the level l-1 basis from the class-r population
 * constraint, as in setup1.  With ntag = 1 and q = (0..d-1) the rows are
 * exactly those of setup1 (the GCE rows are its subtractive rows), so
 * --bmax 1 reproduces the original gmom. */
static void set_q_from_z(mpq_t q, mpz_t z) { mpq_set_z(q, z); }
static void set_q_neg_z(mpq_t q, mpz_t z) { mpq_set_z(q, z); mpq_neg(q, q); }

LSB* setupb(mpz_t** L, int* q, int d, int ntag, const int* clevel, int r, int l, int* N, mpz_t* Z, int aonly)
{
	int i, j, s, row, t, c;
	LSB* ls = (LSB*) calloc(1, sizeof(LSB));

	int szIk = nck(d + l - 1, l);           /* level l   */
	int szI  = nck(d + (l-1) - 1, l-1);     /* level l-1 */
	int** Ik = sortbynnzpos(multichoose(d, l),   szIk, d);
	int** I  = sortbynnzpos(multichoose(d, l-1), szI,  d);
	int colA = szIk * r, colD = szI * r;
	ls->szIk = szIk; ls->szI = szI; ls->colA = colA; ls->colD = colD;

	/* child combination tables: the child drops local queue d-ntag+c and
	 * carries clevel[c] = l or l+1; the GCE term is at level l there */
	int szIi = (d >= 2) ? nck((d-1) + l - 1, l) : 0;
	int** Ii = (d >= 2) ? sortbynnzpos(multichoose(d-1, l), szIi, d-1) : NULL;
	int* tmp = (int*) malloc((d > 0 ? d : 1) * sizeof(int));

	int nrows = 0;
	for (c = 0; c < ntag; c++) { int jj = d - ntag + c; for (i = 0; i < szIk; i++) if (Ik[i][jj] == 0) nrows++; }
	for (i = 0; i < szIk; i++) for (j = 0; j < d; j++) if (Ik[i][j] > 0) nrows++;
	nrows += szI * (r - 1);
	ls->nrows = nrows;

	ls->A = mpq_matzeros(nrows, colA);
	if (!aonly) {
		ls->B = mpq_matzeros(nrows, colA);
		ls->D = mpq_matzeros(nrows, colD);
		ls->E = mpq_matzeros(colD, colA);
		ls->F = mpq_matzeros(colD, colD);
		ls->cchild = (int*) malloc(nrows * sizeof(int));
		ls->cslot  = (int*) malloc(nrows * sizeof(int));
		ls->ccol   = (int*) malloc(nrows * sizeof(int));
		for (i = 0; i < nrows; i++) ls->cchild[i] = -1;
	}

	row = 0;
	/* GCE: queue q[jj] removed, child c */
	for (c = 0; c < ntag; c++) {
		int jj = d - ntag + c, k = q[jj];
		for (i = 0; i < szIk; i++) {
			if (Ik[i][jj] != 0) continue;
			mpq_set_si(ls->A[row][i*r + 0], 1, 1);
			for (s = 1; s < r; s++) set_q_neg_z(ls->A[row][i*r + s], L[k][s-1]);
			if (!aonly) {
				set_q_from_z(ls->B[row][i*r + 0], L[k][r-1]);
				int a, b = 0;
				for (a = 0; a < d; a++) if (a != jj) tmp[b++] = Ik[i][a];
				t = (int) int_matmatchrow(Ii, szIi, d-1, tmp);
				ls->cchild[row] = c;
				ls->cslot[row]  = (clevel[c] == l) ? 0 : 1;   /* 0: vl, 1: vl1 */
				ls->ccol[row]   = t*r + 0;
			}
			row++;
		}
	}

	/* CE: a replica of queue q[j] removed */
	for (i = 0; i < szIk; i++) {
		for (j = 0; j < d; j++) {
			if (Ik[i][j] == 0) continue;
			int k = q[j];
			mpq_set_si(ls->A[row][i*r + 0], 1, 1);
			for (s = 1; s < r; s++) set_q_neg_z(ls->A[row][i*r + s], L[k][s-1]);
			if (!aonly) {
				set_q_from_z(ls->B[row][i*r + 0], L[k][r-1]);
				Ik[i][j]--;
				t = (int) int_matmatchrow(I, szI, d, Ik[i]);
				Ik[i][j]++;
				mpq_set_si(ls->D[row][t*r + 0], 1, 1);
			}
			row++;
		}
	}

	/* PC, classes 1..r-1 */
	for (i = 0; i < szI; i++) {
		for (s = 1; s < r; s++) {
			if (!aonly) {
				mpq_set_si(ls->D[row][i*r + 0], -N[s-1], 1);
				set_q_from_z(ls->D[row][i*r + s], Z[s-1]);
			}
			for (j = 0; j < d; j++) {
				I[i][j]++;
				t = (int) int_matmatchrow(Ik, szIk, d, I[i]);
				I[i][j]--;
				mpz_t v; mpz_init(v);
				mpz_mul_ui(v, L[q[j]][s-1], (unsigned long)(1 + I[i][j]));
				mpz_neg(v, v);
				mpq_set_z(ls->A[row][t*r + s], v);
				mpz_clear(v);
			}
			row++;
		}
	}

	/* PC, class r: level l-1 from the previous population step */
	if (!aonly) {
		row = 0;
		for (i = 0; i < szI; i++) {
			for (s = 0; s < r; s++) {
				set_q_from_z(ls->F[row][i*r + s], Z[r-1]);
				for (j = 0; j < d; j++) {
					I[i][j]++;
					t = (int) int_matmatchrow(Ik, szIk, d, I[i]);
					I[i][j]--;
					mpz_t v; mpz_init(v);
					mpz_mul_ui(v, L[q[j]][r-1], (unsigned long)(1 + I[i][j]));
					mpq_set_z(ls->E[row][t*r + s], v);
					mpz_clear(v);
				}
				row++;
			}
		}
	}

	for (i = 0; i < szIk; i++) free(Ik[i]);
	free(Ik);
	for (i = 0; i < szI; i++) free(I[i]);
	free(I);
	if (Ii) { for (i = 0; i < szIi; i++) free(Ii[i]); free(Ii); }
	free(tmp);
	return ls;
}

static void matfree(mpq_mat_t A, int nr, int nc)
{
	int i, j;
	if (!A) return;
	for (i = 0; i < nr+1; i++) { for (j = 0; j < nc+1; j++) mpq_clear(A[i][j]); free(A[i]); }
	free(A);
}

void lsb_free(LSB* ls)
{
	if (!ls) return;
	matfree(ls->A, ls->nrows, ls->colA);
	matfree(ls->B, ls->nrows, ls->colA);
	matfree(ls->D, ls->nrows, ls->colD);
	matfree(ls->E, ls->colD, ls->colA);
	matfree(ls->F, ls->colD, ls->colD);
	free(ls->cchild); free(ls->cslot); free(ls->ccol);
	free(ls);
}
