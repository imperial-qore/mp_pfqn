#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "promom.h"

/* qdecrease - replica descent for the marginal quantities.
 *
 * Port of mom/mdecrease.c to q_d(P, n).  The descent peels one replica of
 * queue 0 at a time (CE) and fills the class decrements from the population
 * constraints (PC), from level R down to level 0, where the level-0 base
 * component is q(N, n) itself.
 *
 * Two changes with respect to mdecrease:
 *
 *   - the CE is unchanged, since a replica of any station leaves the count
 *     at the reference copy alone;
 *   - the PC gains the n-coupled term, and its station sum runs over the
 *     model with the reference COPY removed, so station M carries the
 *     multiplicity mi_M + i_M - 1.  With Z_s != 0,
 *
 *       q(i, P-1_s, n) = [ N_s q(i, P, n) - n L_Ms q(i, P-1_s, n-1)
 *                          - sum_j m'_j L_js q(i+1_j, P-1_s, n) ] / Z_s,
 *       m'_j = mi_j + i_j - [j = M]
 *
 *     and with Z_s = 0 mom uses the server-count identity
 *     sum_j m_j G(m+1_j, Q) = (|Q| + |m|) G(m, Q) instead.  Its conditioned
 *     form follows by summing the PC over all classes and substituting the
 *     n-recursion sum_s L_Ms q(P-1_s, n-1) = q(P, n):
 *
 *       sum_j m'_j q(m+1_j, Q, n) = (|Q| - n + |m| - 1) q(m, Q, n)
 *
 *     i.e. mom's divisor less the n jobs pinned at the reference copy and
 *     less that copy itself.
 *
 * The n-1 term is at the same level and the same component as the value
 * being computed, so each level is swept with n ascending.
 */
int qdecrease(qnmodel* qnm, mpq_vec_t* gk, mpq_vec_t* g, mpq_vec_t* grv,
              int sumN, mpq_t* dist)
{
	int M = qnm->M, R = qnm->R;
	int t, i, s, l, j, n;
	mpq_t tmp, tmp2;
	mpq_init(tmp); mpq_init(tmp2);

	int Ntot = 0;
	for (s = 0; s < R; s++) Ntot += qnm->N[s];

	/* level-R working set: R+1 components per combination, component R
	 * being the class-R decrement taken from the previous point */
	int cardR = nck(M + R - 1, R);
	mpq_vec_t* Gk = (mpq_vec_t*) calloc(sumN + 1, sizeof(mpq_vec_t));
	for (n = 0; n <= sumN; n++) {
		Gk[n] = mpq_vec(cardR * (R+1), 0, 1);
		for (t = 0; t < cardR; t++) {
			mpq_set(Gk[n][t*(R+1) + R], grv[n][t*R]);
			for (s = 0; s < R; s++)
				mpq_set(Gk[n][t*(R+1) + s], g[n][t*R + s]);
		}
	}
	(void) gk;

	mpq_vec_t* G = NULL;
	for (l = R; l >= 1; l--) {
		int cardl  = nck(M + l - 1, l);
		int cardl1 = nck(M + l - 2, l - 1);
		int** Ik = sortbynnzpos((int**) multichoose(M, l),   cardl,  M);
		int** I  = sortbynnzpos((int**) multichoose(M, l-1), cardl1, M);

		G = (mpq_vec_t*) calloc(sumN + 1, sizeof(mpq_vec_t));
		for (n = 0; n <= sumN; n++) G[n] = mpq_vec(cardl1 * (R+1), 0, 1);

		for (n = 0; n <= sumN; n++) {
			for (i = 0; i < cardl1; i++) {
				/* CE: peel one replica off queue 0 */
				I[i][0]++;
				t = (int) int_matmatchrow(Ik, cardl, M, I[i]);
				I[i][0]--;
				mpq_set(G[n][i*(R+1)], Gk[n][t*(R+1)]);
				for (s = 0; s < R; s++) {
					mpq_set_z(tmp, qnm->L[0][s]);
					mpq_mul(tmp, tmp, Gk[n][t*(R+1) + 1 + s]);
					mpq_sub(G[n][i*(R+1)], G[n][i*(R+1)], tmp);
				}

				/* PC: class decrements */
				for (s = 1; s <= R; s++) {
					if (mpz_cmp_ui(qnm->Z[s-1], 0) == 0) {
						int msum = 0;
						mpq_set_si(G[n][i*(R+1) + s], 0, 1);
						for (j = 0; j < M; j++) {
							I[i][j]++;
							t = (int) int_matmatchrow(Ik, cardl, M, I[i]);
							I[i][j]--;
							mpq_set_si(tmp, qnm->mi[j] + I[i][j] - (j == M-1 ? 1 : 0), 1);
							mpq_mul(tmp, tmp, Gk[n][t*(R+1) + s]);
							mpq_add(G[n][i*(R+1) + s], G[n][i*(R+1) + s], tmp);
							msum += qnm->mi[j];
						}
						for (j = 0; j < M; j++) msum += I[i][j];
						/* mom's (Ntot + msum - 1), less the n pinned jobs and
						 * the reference copy held out of the server count */
						mpq_set_si(tmp, Ntot + msum - 2 - n, 1);
						if (mpq_sgn(tmp) == 0) { mpq_set_si(G[n][i*(R+1) + s], 0, 1); continue; }
						mpq_div(G[n][i*(R+1) + s], G[n][i*(R+1) + s], tmp);
					} else {
						mpq_set_si(G[n][i*(R+1) + s], qnm->N[s-1], 1);
						mpq_set_z(tmp, qnm->Z[s-1]);
						mpq_div(G[n][i*(R+1) + s], G[n][i*(R+1) + s], tmp);
						mpq_mul(G[n][i*(R+1) + s], G[n][i*(R+1) + s], G[n][i*(R+1)]);
						for (j = 0; j < M; j++) {
							I[i][j]++;
							t = (int) int_matmatchrow(Ik, cardl, M, I[i]);
							I[i][j]--;
							mpq_set_z(tmp, qnm->L[j][s-1]);
							mpq_set_si(tmp2, qnm->mi[j] + I[i][j] - (j == M-1 ? 1 : 0), 1);
							mpq_mul(tmp2, tmp2, tmp);
							mpq_set_z(tmp, qnm->Z[s-1]);
							mpq_div(tmp2, tmp2, tmp);
							mpq_mul(tmp2, tmp2, Gk[n][t*(R+1) + s]);
							mpq_sub(G[n][i*(R+1) + s], G[n][i*(R+1) + s], tmp2);
						}
						/* n-coupled term: - n L_Ms q(i, P-1_s, n-1) / Z_s */
						if (n > 0) {
							mpq_set_si(tmp2, n, 1);
							mpq_set_z(tmp, qnm->L[M-1][s-1]);
							mpq_mul(tmp2, tmp2, tmp);
							mpq_set_z(tmp, qnm->Z[s-1]);
							mpq_div(tmp2, tmp2, tmp);
							mpq_mul(tmp2, tmp2, G[n-1][i*(R+1) + s]);
							mpq_sub(G[n][i*(R+1) + s], G[n][i*(R+1) + s], tmp2);
						}
					}
				}
			}
		}
		free(Ik); free(I);

		if (l > 1) {
			for (n = 0; n <= sumN; n++) {
				int len = cardl * (R+1);
				for (t = 0; t < len + 1; t++) mpq_clear(Gk[n][t]);
				free(Gk[n]);
				Gk[n] = G[n];
			}
			free(G);
			G = NULL;
		}
	}

	for (n = 0; n <= sumN; n++) mpq_set(dist[n], G[n][0]);

	for (n = 0; n <= sumN; n++) {
		int len = 1 * (R+1);
		for (t = 0; t < len + 1; t++) mpq_clear(G[n][t]);
		free(G[n]);
	}
	free(G);
	free(Gk);
	mpq_clear(tmp); mpq_clear(tmp2);
	return 0;
}
