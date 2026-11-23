/* gbound.h: a-priori size bounds on the exact normalizing constants that
 * util/singular.h exports (G, G(N-1_r), G^{+k}(N-1_r)), so that a
 * multi-modular solver knows the modulus that suffices before it starts.
 *
 * Denominator.  For integer demands and think times,
 *     G(n) = sum prod_r Z_r^{n0r}/n0r!  prod_k |n_k|!/prod_r n_kr! L_kr^{n_kr}
 * and every station factor is an integer multinomial, so the denominator
 * divides D(n) = prod_{r : Z_r != 0} n_r!  (D = 1 without think times).
 *
 * Value.  G(n) is the coefficient of x^n in
 *     F(x) = prod_r exp(Z_r x_r) prod_k (1 - sum_r L_kr x_r)^(-m_k),
 * whose coefficients are all nonnegative, so for EVERY x > 0 with
 * sum_r L_kr x_r < 1 at each station,  G(n) <= F(x) / x^n.
 * Minimizing over x gives the saddle-point (asymptotic) estimate of
 * log G (SIGMETRICS'17), here used as a bound: an inexact minimum only
 * loosens it.  log F(e^y) - n.y is convex in y; coordinate Newton steps.
 *
 * With 0 <= num <= G_up D and den <= D, zp_ratrecon_bounded(.., D)
 * returns the exact value once M > 2 G_up D^2, with no stability or
 * witness check.  Header-only, used by mommod/engine.c and
 * comommod/engine.c. */
#ifndef GBOUND_H
#define GBOUND_H

#include <math.h>
#include <gmp.h>
#include "util.h"

/* upper bound on log2 G(n) of q with m_k = mi_k (+1 at station kx if
 * kx >= 0); -INFINITY when G(n) = 0 (a negative population, or a class
 * with jobs and nowhere to go) */
static double gb_log2_G(const qnmodel* q, const int* n, int kx)
{
	int M = q->M, R = q->R, k, r, it;
	for (r = 0; r < R; r++) if (n[r] < 0) return -INFINITY;
	double* L = (double*) malloc((size_t)M * R * sizeof(double));
	double* Z = (double*) malloc(R * sizeof(double));
	double* m = (double*) malloc(M * sizeof(double));
	double* y = (double*) malloc(R * sizeof(double));
	double* s = (double*) malloc(M * sizeof(double));
	int*    act = (int*) malloc(R * sizeof(int));
	double res = 0.0;
	long ntot = 0;
	for (r = 0; r < R; r++) { Z[r] = mpz_get_d(q->Z[r]); ntot += n[r]; }
	for (k = 0; k < M; k++) {
		m[k] = (double)q->mi[k] / qn_mden(q) + (k == kx ? 1 : 0);
		for (r = 0; r < R; r++) L[k*R+r] = mpz_get_d(q->L[k][r]);
	}
	/* classes without jobs take x_r = 0: their factors are then 1 */
	for (r = 0; r < R; r++) {
		act[r] = n[r] > 0;
		if (!act[r]) continue;
		double w = Z[r];
		for (k = 0; k < M; k++) w += m[k] * L[k*R+r];
		if (w <= 0.0) { res = -INFINITY; goto out; }
		y[r] = log((double)n[r] / (2.0 * (double)ntot * w));   /* then sum_r L_kr x_r <= 1/2 */
	}
	if (ntot == 0) { res = 0.0; goto out; }
	for (k = 0; k < M; k++) { s[k] = 0; for (r = 0; r < R; r++) if (act[r]) s[k] += L[k*R+r] * exp(y[r]); }
	for (it = 0; it < 200; it++) {
		double moved = 0.0;
		for (r = 0; r < R; r++) {
			if (!act[r]) continue;
			double x = exp(y[r]), a = Z[r], b = 0.0, c;
			for (k = 0; k < M; k++) if (m[k] > 0 && L[k*R+r] > 0) {
				double u = L[k*R+r] / (1.0 - s[k]);
				a += m[k] * u; b += m[k] * u * u;
			}
			double g = x * a - n[r];                 /* d f / d y_r     */
			double h = x * a + x * x * b;            /* d^2 f / d y_r^2 */
			double step = -g / h, f0 = 0.0, f1;
			/* f restricted to y_r, relative to the current point */
			for (c = 1.0; c > 1e-12; c *= 0.5) {
				double xn = exp(y[r] + c * step), dx = xn - x, ok = 1;
				f1 = Z[r] * dx - n[r] * c * step;
				for (k = 0; k < M && ok; k++) if (m[k] > 0 && L[k*R+r] > 0) {
					double sn = s[k] + L[k*R+r] * dx;
					if (sn >= 1.0 - 1e-12) ok = 0;
					else f1 -= m[k] * (log1p(-sn) - log1p(-s[k]));
				}
				if (ok && f1 <= f0) break;
			}
			if (c <= 1e-12) continue;
			{
				double xn = exp(y[r] + c * step), dx = xn - x;
				for (k = 0; k < M; k++) s[k] += L[k*R+r] * dx;
				y[r] += c * step;
				if (fabs(c * step) > moved) moved = fabs(c * step);
			}
		}
		if (moved < 1e-9) break;
	}
	/* f(y) = sum Z_r x_r - sum m_k log(1 - s_k) - sum n_r y_r, each term
	 * recomputed from y; the margin covers the double rounding of a sum of
	 * terms of this size with room to spare */
	{
		double f = 0.0, mag = 0.0, t;
		for (k = 0; k < M; k++) {
			double sk = 0.0;
			for (r = 0; r < R; r++) if (act[r]) sk += L[k*R+r] * exp(y[r]);
			if (m[k] == 0 || sk == 0.0) continue;
			t = -m[k] * log1p(-sk); f += t; mag += fabs(t);
		}
		for (r = 0; r < R; r++) if (act[r]) {
			t = Z[r] * exp(y[r]) - n[r] * y[r]; f += t;
			mag += fabs(Z[r] * exp(y[r])) + fabs(n[r] * y[r]);
		}
		res = f / log(2.0) + 1e-9 * mag + 8.0;
	}
out:
	free(L); free(Z); free(m); free(y); free(s); free(act);
	return res;
}

/* D(n) = prod_{r : Z_r != 0} n_r!, a multiple of the denominator of G(n) */
static void gb_den(mpz_t D, const qnmodel* q, const int* n)
{
	int r; mpz_t f; mpz_init(f);
	mpz_set_ui(D, 1);
	for (r = 0; r < q->R; r++)
		if (mpz_sgn(q->Z[r]) && n[r] > 0) { mpz_fac_ui(f, (unsigned long)n[r]); mpz_mul(D, D, f); }
	mpz_clear(f);
}

/* Bounds for the 1 + R + M R exported constants of q (util/singular.h
 * layout): D[i] bounds the denominator of quantity i (all initialized
 * here, the caller clears them).
 * Returns the modulus size in bits that makes every bounded reconstruction
 * exact, or -1 if a bound is not finite. */
static long gb_modbits(const qnmodel* q, mpz_t* D)
{
	int M = q->M, R = q->R, r, i;
	int* n = (int*) malloc(R * sizeof(int));
	double need = 0.0;
	for (i = 0; i < 1 + R + M * R; i++) mpz_init(D[i]);   /* all of them, even on failure */
	for (i = 0; i < 1 + R + M * R; i++) {
		int s = -1, kx = -1;
		if (i == 0) ;
		else if (i <= R) s = i - 1;
		else { s = (i - 1 - R) % R; kx = (i - 1 - R) / R; }
		for (r = 0; r < R; r++) n[r] = q->N[r] - (r == s);
		gb_den(D[i], q, n);
		double lg = gb_log2_G(q, n, kx);
		if (isinf(lg) && lg < 0) continue;     /* G = 0: reconstructs from any modulus */
		if (!isfinite(lg)) { free(n); return -1; }
		/* M > 2 G_up D^2 */
		double b = 1.0 + (lg > 0 ? lg : 0.0) + 2.0 * (double)mpz_sizeinbase(D[i], 2);
		if (b > need) need = b;
	}
	free(n);
	return (long)ceil(need) + 2;
}

#endif
