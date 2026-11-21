#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

/* Print G(N), throughputs X_r and queue lengths Q_kr from levels 0 and 1
 * of the full gmom basis at N, in the layout
 *   G[0] = G(N), G[s] = G(N-1_s),
 *   Gk[(k-1)*(R+1)+s] = G(M+1_k, N-1_s)   (s=0 is N).
 * main.c carries the lower levels through the recursion (util/lvtab.h);
 * the replica descent that used to rebuild them here used, for Z_s=0, an
 * identity that fails when other classes have positive think times. */
void gmom_measures(qnmodel* qn, mpq_vec_t G, mpq_vec_t Gk,
                   int out_e, int out_g, int out_l, int out_t, int out_q, int* cperm, mpz_t scale_factor)
{
	int M = qn->M, R = qn->R;
	/* cperm[j] = original class at permuted position j; inv[c] = permuted
	 * position of original class c.  X/Q are printed in ORIGINAL order. */
	int inv[64], _c; for (_c = 0; _c < R; _c++) inv[_c] = _c;
	if (cperm) for (_c = 0; _c < R; _c++) inv[cperm[_c]] = _c;
	mpz_t** L = qn->L; int* N = qn->N; int* mi = qn->mi;
	int i, s;
	mpq_t tmp, tmp2; mpq_init(tmp); mpq_init(tmp2);

	int Ntot = 0; for (s = 0; s < R; s++) Ntot += N[s];

	/* Undo the perturbation scaling (mom/perfindices convention):
	 *   G      is degree Ntot in the demands  -> divide by scale^Ntot;
	 *   X_r    = G(N-1_r)/G(N)                 -> multiply by scale;
	 *   Q_kr   the two scale powers cancel     -> no correction, but a
	 *          demand that was 0 before perturbing must give Q=0. */
	int scaled = (mpz_cmp_ui(scale_factor, 1) > 0);
	mpq_t G0s; mpq_init(G0s); mpq_set(G0s, G[0]);
	mpq_t qscale; mpq_init(qscale); mpq_set_z(qscale, scale_factor);
	if (scaled) {
		mpz_t sp; mpz_init(sp); mpz_pow_ui(sp, scale_factor, Ntot);
		mpq_t d; mpq_init(d); mpq_set_z(d, sp);
		mpq_div(G0s, G0s, d);
		mpq_clear(d); mpz_clear(sp);
	}

	if (out_l) {
		printf("%.15e\n", mpq_log(G0s)); /* no overflow for G > 1e308 */
	} else if (out_e) {
		mpq_canonicalize(G0s);
		mpz_t num, den; mpz_init(num); mpz_init(den);
		mpq_get_num(num, G0s); mpq_get_den(den, G0s);
		gmp_printf("%Zd\n%Zd\n", num, den);
		mpz_clear(num); mpz_clear(den);
	} else if (out_g) {
		printf("%.15e\n", mpq_get_d(G0s));
	} else if (out_t) {
		for (_c = 0; _c < R; _c++) { mpq_div(tmp, G[inv[_c]+1], G[0]); if (scaled) mpq_mul(tmp, tmp, qscale); printf("%.15e\n", mpq_get_d(tmp)); }
	} else if (out_q) {
		for (i = 1; i <= M; i++) {
			for (_c = 0; _c < R; _c++) {
				s = inv[_c] + 1;
				mpz_t oL; mpz_init(oL); if (scaled) mpz_tdiv_q(oL, L[i-1][s-1], scale_factor); else mpz_set(oL, L[i-1][s-1]);
				if (mpz_sgn(oL) == 0) mpq_set_si(tmp, 0, 1);
				else { mpq_set_z(tmp, L[i-1][s-1]); mpq_mul(tmp, tmp, Gk[(i-1)*(R+1)+s]); mpq_div(tmp, tmp, G[0]); mpq_set_si(tmp2, mi[i-1], 1); mpq_mul(tmp, tmp, tmp2); }
				mpz_clear(oL);
				printf("%.15e", mpq_get_d(tmp));
				if (_c < R-1) printf(" ");
			}
			printf("\n");
		}
	} else {
		printf("========== gmom (generalized MoM, b=1) ==========\n");
		if (scaled) printf("(perturbed approximation)\n");
		printf("G = %.15e\n", mpq_get_d(G0s));
		printf("\nX (throughputs):\n");
		for (_c = 0; _c < R; _c++) { mpq_div(tmp, G[inv[_c]+1], G[0]); if (scaled) mpq_mul(tmp, tmp, qscale); printf("X[%d] = %.15e\n", _c+1, mpq_get_d(tmp)); }
		printf("\nQ (mean queue lengths):\n");
		for (i = 1; i <= M; i++) {
			printf("Q[%d] =", i);
			for (_c = 0; _c < R; _c++) {
				s = inv[_c] + 1;
				mpz_t oL; mpz_init(oL); if (scaled) mpz_tdiv_q(oL, L[i-1][s-1], scale_factor); else mpz_set(oL, L[i-1][s-1]);
				if (mpz_sgn(oL) == 0) mpq_set_si(tmp, 0, 1);
				else { mpq_set_z(tmp, L[i-1][s-1]); mpq_mul(tmp, tmp, Gk[(i-1)*(R+1)+s]); mpq_div(tmp, tmp, G[0]); mpq_set_si(tmp2, mi[i-1], 1); mpq_mul(tmp, tmp, tmp2); }
				mpz_clear(oL);
				printf("\t%.15e", mpq_get_d(tmp));
			}
			printf("\n");
		}
		printf("=========================================\n");
	}

	mpq_clear(tmp); mpq_clear(tmp2); mpq_clear(G0s); mpq_clear(qscale);
}
