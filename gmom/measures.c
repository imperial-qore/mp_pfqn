#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

/* Recover plain G(N), throughputs X_r and queue lengths Q_kr from the
 * gmom moment basis, by the mdecrease replica-descent of the original MoM
 * (mom/mdecrease.c) plus perfindices.
 *
 * gmom's top basis vlM is the level-(R-1) basis of the full M-queue model:
 * vlM[t*R+s] = G(M+delta_t, N-1_s) for delta_t = multichoose(M,R-1)[t] and
 * s=0..R-1 (s=0 is N).  grM is the same basis one population step back,
 * so grM[t*R+0] = G(M+delta_t, N-1_R) supplies the missing class-R
 * decrement.  The descent peels replicas off queue 0 via the convolution
 * expression and fills the class decrements via the population
 * constraints, exactly as the original MoM, but starts one level lower
 * (R-1 instead of R) because the generalized basis is more compact.
 */
void gmom_measures(qnmodel* qn, mpq_vec_t vlM, mpq_vec_t grM,
                   int out_e, int out_g, int out_l, int out_t, int out_q)
{
	int M = qn->M, R = qn->R;
	mpz_t** L = qn->L; mpz_t* Z = qn->Z; int* N = qn->N; int* mi = qn->mi;
	int t, i, s, l, j;
	mpq_t tmp, tmp2; mpq_init(tmp); mpq_init(tmp2);

	int Ntot = 0; for (s = 0; s < R; s++) Ntot += N[s];

	/* level-(R-1) Gk with R+1 components */
	int L1 = R - 1;
	int nc = nck(M + L1 - 1, L1);
	mpq_vec_t Gk = mpq_vec(nc * (R+1), 0, 1);
	for (t = 0; t < nc; t++) {
		mpq_set(Gk[t*(R+1) + R], grM[t*R]);
		for (s = 0; s < R; s++) mpq_set(Gk[t*(R+1) + s], vlM[t*R + s]);
	}

	mpq_vec_t G = NULL;
	for (l = L1; l >= 1; l--) {
		int ncl  = nck(M + l - 1, l);      /* combos at level l   */
		int ncl1 = nck(M + l - 2, l - 1);  /* combos at level l-1 */
		G = mpq_vec(ncl1 * (R+1), 0, 1);
		int** Ik = sortbynnzpos(multichoose(M, l),   ncl,  M);
		int** I  = sortbynnzpos(multichoose(M, l-1), ncl1, M);
		for (i = 0; i < ncl1; i++) {
			/* CE: peel a replica off queue 0 */
			I[i][0]++;
			t = (int) int_matmatchrow(Ik, ncl, M, I[i]);
			I[i][0]--;
			mpq_set(G[i*(R+1)], Gk[t*(R+1)]);
			for (s = 0; s < R; s++) {
				mpq_set_z(tmp, L[0][s]);
				mpq_mul(tmp, tmp, Gk[t*(R+1) + 1 + s]);
				mpq_sub(G[i*(R+1)], G[i*(R+1)], tmp);
			}
			/* PC: class decrements */
			for (s = 1; s <= R; s++) {
				if (mpz_cmp_ui(Z[s-1], 0) == 0) {
					mpq_set_si(G[i*(R+1) + s], 0, 1);
					int msum = 0;
					for (j = 0; j < M; j++) {
						I[i][j]++;
						t = (int) int_matmatchrow(Ik, ncl, M, I[i]);
						I[i][j]--;
						mpq_set_si(tmp, mi[j] + I[i][j], 1);
						mpq_mul(tmp, tmp, Gk[t*(R+1) + s]);
						mpq_add(G[i*(R+1) + s], G[i*(R+1) + s], tmp);
						msum += mi[j];
					}
					for (j = 0; j < M; j++) msum += I[i][j];
					mpq_set_si(tmp, Ntot + msum - 1, 1);
					mpq_div(G[i*(R+1) + s], G[i*(R+1) + s], tmp);
				} else {
					mpq_set_si(G[i*(R+1) + s], N[s-1], 1);
					mpq_set_z(tmp, Z[s-1]);
					mpq_div(G[i*(R+1) + s], G[i*(R+1) + s], tmp);
					mpq_mul(G[i*(R+1) + s], G[i*(R+1) + s], G[i*(R+1)]);
					for (j = 0; j < M; j++) {
						I[i][j]++;
						t = (int) int_matmatchrow(Ik, ncl, M, I[i]);
						I[i][j]--;
						mpq_set_z(tmp, L[j][s-1]);
						mpq_set_si(tmp2, mi[j] + I[i][j], 1);
						mpq_mul(tmp2, tmp2, tmp);
						mpq_set_z(tmp, Z[s-1]);
						mpq_div(tmp2, tmp2, tmp);
						mpq_mul(tmp2, tmp2, Gk[t*(R+1) + s]);
						mpq_sub(G[i*(R+1) + s], G[i*(R+1) + s], tmp2);
					}
				}
			}
		}
		free(Ik); free(I);
		if (l > 1) {
			mpq_vecdup(Gk, G, ncl1 * (R+1));
			for (i = 0; i < ncl1*(R+1)+1; i++) mpq_clear(G[i]);
			free(G);
		}
	}
	/* G is level-0: G[0]=G(N), G[s]=G(N-1_s).  Gk is level-1 (single
	 * replica per queue): Gk[(k-1)*(R+1)+r] = G(M+1_k, N-1_r). */

	if (out_l) {
		mpf_t f; mpf_init(f); mpf_set_q(f, G[0]);
		printf("%.15e\n", log(mpf_get_d(f))); mpf_clear(f);
	} else if (out_e) {
		mpq_canonicalize(G[0]);
		mpz_t num, den; mpz_init(num); mpz_init(den);
		mpq_get_num(num, G[0]); mpq_get_den(den, G[0]);
		gmp_printf("%Zd\n%Zd\n", num, den);
		mpz_clear(num); mpz_clear(den);
	} else if (out_g) {
		printf("%.15e\n", mpq_get_d(G[0]));
	} else if (out_t) {
		for (s = 1; s <= R; s++) { mpq_div(tmp, G[s], G[0]); printf("%.15e\n", mpq_get_d(tmp)); }
	} else if (out_q) {
		for (i = 1; i <= M; i++) {
			for (s = 1; s <= R; s++) {
				mpq_set_z(tmp, L[i-1][s-1]);
				mpq_mul(tmp, tmp, Gk[(i-1)*(R+1) + s]);
				mpq_div(tmp, tmp, G[0]);
				mpq_set_si(tmp2, mi[i-1], 1); mpq_mul(tmp, tmp, tmp2);
				printf("%.15e", mpq_get_d(tmp));
				if (s < R) printf(" ");
			}
			printf("\n");
		}
	} else {
		printf("========== gmom (generalized MoM, b=1) ==========\n");
		printf("G = %.15e\n", mpq_get_d(G[0]));
		printf("\nX (throughputs):\n");
		for (s = 1; s <= R; s++) { mpq_div(tmp, G[s], G[0]); printf("X[%d] = %.15e\n", s, mpq_get_d(tmp)); }
		printf("\nQ (mean queue lengths):\n");
		for (i = 1; i <= M; i++) {
			printf("Q[%d] =", i);
			for (s = 1; s <= R; s++) {
				mpq_set_z(tmp, L[i-1][s-1]); mpq_mul(tmp, tmp, Gk[(i-1)*(R+1)+s]); mpq_div(tmp, tmp, G[0]);
				mpq_set_si(tmp2, mi[i-1], 1); mpq_mul(tmp, tmp, tmp2);
				printf("\t%.15e", mpq_get_d(tmp));
			}
			printf("\n");
		}
		printf("=========================================\n");
	}

	mpq_clear(tmp); mpq_clear(tmp2);
}
