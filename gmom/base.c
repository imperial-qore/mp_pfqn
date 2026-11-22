#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

/* Exact normalizing constant of a closed product-form model composed of
 * nq load-independent single-server queues (demand rows dem[q][0..r-1])
 * plus one infinite-server delay station with think times Z[0..r-1], at
 * population pop[0..r-1].  Computed by Buzen convolution over the whole
 * population box, so it is exact for Z=0 and Z!=0 alike.  This is the
 * gmva() the reference MATLAB calls for the base levels of the recursion.
 *
 * Cost is O(box * (nq+1) * r) with box = prod_c (pop[c]+1); used only for
 * the single-queue prefix and the single-class initialisation, where the
 * box is small.
 */

/* mixed-radix helpers over the box [0..pop[0]] x ... x [0..pop[r-1]] */
static int box_size(int* pop, int r)
{
	int t, s = 1;
	for (t = 0; t < r; t++) s *= (pop[t] + 1);
	return s;
}

void gmva(mpq_t G, mpz_t** dem, int nq, int* pop, mpz_t* Z, int r)
{
	int bs = box_size(pop, r);
	int i, c, q, idx;
	int* radix = (int*) malloc(r * sizeof(int));
	int* n = (int*) malloc(r * sizeof(int));
	mpq_vec_t g = mpq_vec(bs, 0, 1);
	mpq_t term, dq;
	mpq_init(term); mpq_init(dq);

	radix[0] = 1;
	for (c = 1; c < r; c++) radix[c] = radix[c-1] * (pop[c-1] + 1);

	/* Station 0: the IS delay.  g[n] = prod_c Z[c]^{n_c} / n_c!.
	 * With Z all zero this collapses to g[0]=1, rest 0. */
	for (idx = 0; idx < bs; idx++) {
		int rem = idx;
		mpq_set_si(g[idx], 1, 1);
		for (c = 0; c < r; c++) { n[c] = rem % (pop[c]+1); rem /= (pop[c]+1); }
		for (c = 0; c < r; c++) {
			if (n[c] == 0) continue;
			if (mpz_sgn(Z[c]) == 0) { mpq_set_si(g[idx], 0, 1); break; }
			/* multiply by Z[c]^{n_c} / n_c! */
			mpq_set_z(dq, Z[c]);
			int e;
			for (e = 0; e < n[c]; e++) {
				mpq_mul(g[idx], g[idx], dq);
				mpq_t inv; mpq_init(inv); mpq_set_si(inv, 1, e+1);
				mpq_mul(g[idx], g[idx], inv); mpq_clear(inv);
			}
		}
	}

	/* Convolve in each LI queue via Buzen: for n ascending,
	 * g_new(n) = g_old(n) + sum_c dem_q[c] * g_new(n - e_c). */
	for (q = 0; q < nq; q++) {
		for (idx = 0; idx < bs; idx++) {
			int rem = idx;
			for (c = 0; c < r; c++) { n[c] = rem % (pop[c]+1); rem /= (pop[c]+1); }
			for (c = 0; c < r; c++) {
				if (n[c] == 0) continue;
				if (mpz_sgn(dem[q][c]) == 0) continue;
				mpq_set_z(dq, dem[q][c]);
				mpq_mul(term, dq, g[idx - radix[c]]);
				mpq_add(g[idx], g[idx], term);
			}
		}
	}

	mpq_set(G, g[bs-1]);

	for (i = 0; i < bs+1; i++) mpq_clear(g[i]);
	free(g);
	mpq_clear(term); mpq_clear(dq);
	free(radix); free(n);
}

/* Convenience wrapper used by the m=1 prefix: `mult` identical copies of a
 * single queue with demand row d, plus delays Z, at population pop. */
void gmom_base(mpq_t G, mpz_t* d, int mult, int* pop, mpz_t* Z, int r)
{
	int q;
	mpz_t** dem = (mpz_t**) malloc(mult * sizeof(mpz_t*));
	for (q = 0; q < mult; q++) dem[q] = d; /* alias: same demand row */
	gmva(G, dem, mult, pop, Z, r);
	free(dem);
}
