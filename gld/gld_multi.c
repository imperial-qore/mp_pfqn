#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gld.h"

/*
 * Multi-class load-dependent normalizing constant via station-by-station
 * convolution (generalized Buzen for load-dependent stations).
 *
 * Single-station balance factor for station m holding class-vector k
 * (total j = sum_r k_r):
 *   Y_m(k) = j! / prod_r(k_r!) * prod_r L[m][r]^{k_r} / prod_{i=1}^{j} mu[m][i-1]
 *
 * Convolution over stations (n ranges over 0 <= n_r <= N_r):
 *   G_0(n) = 1 if n == 0 else 0
 *   G_m(n) = sum_{0 <= k <= n} Y_m(k) * G_{m-1}(n - k)
 *   G(N)   = G_M(N)
 *
 * Think time Z is modelled as one extra delay (infinite-server) station with
 * demand Z_r and rate mu[k] = k+1, whose factor reduces to
 *   Y_Z(k) = prod_r Z_r^{k_r} / prod_r(k_r!).
 * It is convolved in after the M queueing stations when any Z_r > 0.
 *
 * A load-independent station (mu[m][k] == 1 for all k, i.e. a single server)
 * obeys the classical Buzen recursion
 *   G_m(n) = G_{m-1}(n) + sum_r L[m][r] * G_m(n - e_r),
 * an in-place O(P*R) update (P = prod_r (N_r+1)); it is used whenever it
 * applies. Only genuinely load-dependent stations use the full convolution.
 *
 * This replaces the previous recursion, which had no memoization and ran in
 * time exponential in the total population (it hung on models with large Nt).
 */

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */
void gld_auto_mu(qnmodel *qn)
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

void gld_multi(mpq_t G, mpz_t **L, int *N, mpq_t **mu, mpz_t *Z, int M, int R, int Nt)
{
	int m, r, j;
	long int i;

	/* Base case: all populations zero */
	if (Nt == 0) {
		mpq_set_ui(G, 1, 1);
		return;
	}
	/* Base case: no stations and no think time, positive population */
	int hasZ = 0;
	for (r = 0; r < R; r++)
		if (mpz_sgn(Z[r]) > 0) {
			hasZ = 1;
			break;
		}
	if (M == 0 && !hasZ) {
		mpq_set_ui(G, 0, 1);
		return;
	}

	/* Mixed-radix strides over the box 0 <= n_r <= N_r.
	 * P = prod_r (N_r + 1) is the number of population vectors. */
	long int P = 1;
	long int *stride = (long int *)malloc(R * sizeof(long int));
	for (r = 0; r < R; r++) {
		stride[r] = P;
		P *= (long int)(N[r] + 1);
	}

	/* Does any queueing station need the full load-dependent convolution? */
	int has_ld = 0;
	for (m = 0; m < M && !has_ld; m++)
		for (j = 0; j < Nt; j++)
			if (mpq_cmp_ui(mu[m][j], 1, 1) != 0) {
				has_ld = 1;
				break;
			}
	int need_scratch = has_ld || hasZ;

	/* Running normalizing constant table G_m; starts at G_0 = delta_0 */
	mpq_t *Gv = (mpq_t *)malloc(P * sizeof(mpq_t));
	for (i = 0; i < P; i++)
		mpq_init(Gv[i]);
	mpq_set_ui(Gv[0], 1, 1);

	/* Scratch buffers only needed for full convolutions */
	mpq_t *Gs = NULL, *Y = NULL;
	if (need_scratch) {
		Gs = (mpq_t *)malloc(P * sizeof(mpq_t));
		Y = (mpq_t *)malloc(P * sizeof(mpq_t));
		for (i = 0; i < P; i++) {
			mpq_init(Gs[i]);
			mpq_init(Y[i]);
		}
	}

	int *k = (int *)malloc(R * sizeof(int));
	int *nn = (int *)malloc(R * sizeof(int));
	mpz_t num, den, powv;
	mpz_init(num);
	mpz_init(den);
	mpz_init(powv);
	mpq_t Yq, tmp, Lq;
	mpq_init(Yq);
	mpq_init(tmp);
	mpq_init(Lq);

	/* Total number of stations to convolve, including the optional delay */
	int Mtot = M + (hasZ ? 1 : 0);

	for (m = 0; m < Mtot; m++) {
		int is_delay = (m == M); /* the appended delay station */

		if (!is_delay) {
			/* Is this queueing station load-independent (single server)? */
			int li = 1;
			for (j = 0; j < Nt; j++)
				if (mpq_cmp_ui(mu[m][j], 1, 1) != 0) {
					li = 0;
					break;
				}

			if (li) {
				/* In-place Buzen: G[n] += sum_r L[m][r] * G[n - e_r], ascending */
				for (r = 0; r < R; r++)
					nn[r] = 0;
				for (i = 1; i < P; i++) {
					for (r = 0; r < R; r++) {
						if (nn[r] < N[r]) {
							nn[r]++;
							break;
						}
						nn[r] = 0;
					}
					for (r = 0; r < R; r++) {
						if (nn[r] > 0) {
							mpq_set_z(Lq, L[m][r]);
							mpq_mul(tmp, Lq, Gv[i - stride[r]]);
							mpq_add(Gv[i], Gv[i], tmp);
						}
					}
				}
				continue;
			}
		}

		/* Full convolution: build the station factor table Y[.] */
		mpq_t *muprod = NULL;
		if (!is_delay) {
			/* muprod[j] = prod_{i=0}^{j-1} mu[m][i] */
			muprod = (mpq_t *)malloc((Nt + 1) * sizeof(mpq_t));
			mpq_init(muprod[0]);
			mpq_set_ui(muprod[0], 1, 1);
			for (j = 1; j <= Nt; j++) {
				mpq_init(muprod[j]);
				mpq_mul(muprod[j], muprod[j - 1], mu[m][j - 1]);
			}
		}

		for (r = 0; r < R; r++)
			k[r] = 0;
		for (i = 0; i < P; i++) {
			j = 0;
			for (r = 0; r < R; r++)
				j += k[r];

			/* num = j! * prod_r dem[r]^{k_r} for a queueing station;
			 * num = prod_r Z[r]^{k_r} for the delay station (no j! factor,
			 * since its mu-product j! cancels the multinomial's j!). */
			int iszero = 0;
			if (is_delay)
				mpz_set_ui(num, 1);
			else
				mpz_fac_ui(num, j);
			for (r = 0; r < R; r++) {
				if (k[r] > 0) {
					mpz_ptr dem = is_delay ? Z[r] : L[m][r];
					if (mpz_sgn(dem) == 0) {
						iszero = 1;
						break;
					}
					mpz_pow_ui(powv, dem, k[r]);
					mpz_mul(num, num, powv);
				}
			}
			if (iszero) {
				mpq_set_ui(Y[i], 0, 1);
			} else {
				/* den = prod_r k_r! */
				mpz_set_ui(den, 1);
				for (r = 0; r < R; r++) {
					mpz_fac_ui(powv, k[r]);
					mpz_mul(den, den, powv);
				}
				mpq_set_z(Yq, num);
				mpq_set_z(tmp, den);
				mpq_div(Yq, Yq, tmp);
				if (!is_delay)
					mpq_div(Yq, Yq, muprod[j]);
				mpq_set(Y[i], Yq);
			}

			for (r = 0; r < R; r++) {
				if (k[r] < N[r]) {
					k[r]++;
					break;
				}
				k[r] = 0;
			}
		}

		/* Gs(n) = sum_{0 <= kk <= n} Y[kk] * Gv[n - kk] */
		for (r = 0; r < R; r++)
			nn[r] = 0;
		for (i = 0; i < P; i++) {
			mpq_set_ui(Gs[i], 0, 1);
			for (r = 0; r < R; r++)
				k[r] = 0;
			while (1) {
				long int ik = 0, idiff = 0;
				for (r = 0; r < R; r++) {
					ik += (long int)k[r] * stride[r];
					idiff += (long int)(nn[r] - k[r]) * stride[r];
				}
				if (mpq_sgn(Y[ik]) != 0 && mpq_sgn(Gv[idiff]) != 0) {
					mpq_mul(tmp, Y[ik], Gv[idiff]);
					mpq_add(Gs[i], Gs[i], tmp);
				}
				int carry = 1;
				for (r = 0; r < R; r++) {
					if (k[r] < nn[r]) {
						k[r]++;
						carry = 0;
						break;
					}
					k[r] = 0;
				}
				if (carry)
					break;
			}
			for (r = 0; r < R; r++) {
				if (nn[r] < N[r]) {
					nn[r]++;
					break;
				}
				nn[r] = 0;
			}
		}

		for (i = 0; i < P; i++)
			mpq_set(Gv[i], Gs[i]);

		if (!is_delay) {
			for (j = 0; j <= Nt; j++)
				mpq_clear(muprod[j]);
			free(muprod);
		}
	}

	/* idx(N) = sum_r N_r * stride_r = P - 1 */
	mpq_set(G, Gv[P - 1]);

	/* Cleanup */
	mpz_clear(num);
	mpz_clear(den);
	mpz_clear(powv);
	mpq_clear(Yq);
	mpq_clear(tmp);
	mpq_clear(Lq);
	for (i = 0; i < P; i++)
		mpq_clear(Gv[i]);
	free(Gv);
	if (need_scratch) {
		for (i = 0; i < P; i++) {
			mpq_clear(Gs[i]);
			mpq_clear(Y[i]);
		}
		free(Gs);
		free(Y);
	}
	free(stride);
	free(k);
	free(nn);
}
