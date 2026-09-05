#include <stdio.h>
#include <stdlib.h>
#include "clw.h"
#include <math.h>
#include <complex.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* read-only context shared across the recursion */
typedef struct {
	int qd, p;
	const int *N;         /* p */
	const int *l;         /* p */
	const double *r;      /* p : contour radii */
	const double *arho0;  /* p : alpha_j rho_{j0} */
	const double *rhoS;   /* qd*p : alpha_j rho_{ji}, row i chain j */
	const double *m;      /* qd : queue multiplicities */
} clw_ctx;

/* scaled generating function Gbar(w) at the p-vector w
 *   Gbar(w) = exp( sum_j arho0_j (w_j - 1) )
 *             / prod_i (1 - sum_j rhoS_{ji} w_j)^{m_i}
 * (principal complex log, branch choices cancel through exp). */
static double complex clw_gbar(const clw_ctx *c, const double complex *w)
{
	double complex expo = 0.0;
	for (int j = 0; j < c->p; j++)
		expo += c->arho0[j] * (w[j] - 1.0);
	double complex logden = 0.0;
	for (int i = 0; i < c->qd; i++) {
		double complex a = 0.0;
		for (int j = 0; j < c->p; j++)
			a += c->rhoS[i * c->p + j] * w[j];
		logden += c->m[i] * clog(1.0 - a);
	}
	return cexp(expo - logden);
}

/* one-dimensional lattice-Poisson inversion (eq. 2.3), scaled: extracts the
 * coefficient of w_j^{K_j}, recursing on inner chains. w[0..j-1] are fixed. */
static double complex clw_invert(const clw_ctx *c, int j, double complex *w)
{
	int Kj = c->N[j], lj = c->l[j];
	double rj = c->r[j];
	double complex acc = 0.0;
	for (int k1 = 0; k1 < lj; k1++) {
		double complex ph = cexp(-I * M_PI * k1 / lj);
		double complex inner = 0.0;
		for (int k = -Kj; k < Kj; k++) {
			double sign = (k % 2 == 0) ? 1.0 : -1.0;
			double theta = M_PI * (k1 + (double)lj * k) / ((double)lj * Kj);
			w[j] = rj * cexp(I * theta);
			double complex v = (j == c->p - 1) ? clw_gbar(c, w)
			                                   : clw_invert(c, j + 1, w);
			inner += sign * v;
		}
		acc += ph * inner;
	}
	double complex val = acc / (2.0 * lj * Kj * pow(rj, Kj));
	if (j == 0)
		val = creal(val);
	return val;
}

void pfqn_clw(int qd, int p0, const double *L0, const int *N0, const double *Z0,
              const double *m0, const int *l0, const double *gam0,
              double *Gout, double *lGout)
{
	/* trivial populations */
	int anyNeg = 0, sumN = 0;
	for (int j = 0; j < p0; j++) {
		if (N0[j] < 0) anyNeg = 1;
		if (N0[j] > 0) sumN += N0[j];
	}
	if (anyNeg) { *Gout = 0.0; *lGout = -INFINITY; return; }
	if (sumN == 0) { *Gout = 1.0; *lGout = 0.0; return; }

	/* drop zero-population chains: the coefficient of z_j^0 is the pgf
	 * restricted to z_j = 0, which removes chain j exactly. */
	int p = 0;
	int *keep = (int *)malloc((size_t)p0 * sizeof(int));
	for (int j = 0; j < p0; j++) {
		keep[j] = (N0[j] > 0);
		if (keep[j]) p++;
	}
	int *N = (int *)malloc((size_t)p * sizeof(int));
	int *l = (int *)malloc((size_t)p * sizeof(int));
	double *gam = (double *)malloc((size_t)p * sizeof(double));
	double *Z = (double *)malloc((size_t)p * sizeof(double));
	double *L = (double *)malloc((size_t)qd * p * sizeof(double));
	for (int j = 0, jj = 0; j < p0; j++) {
		if (!keep[j]) continue;
		N[jj] = N0[j];
		Z[jj] = Z0 ? Z0[j] : 0.0;
		/* CLW default lattice/aliasing parameters (Section 2.2, p. 962) */
		if (l0) l[jj] = l0[j];
		else    l[jj] = (jj == 0) ? 1 : (jj <= 2 ? 2 : 3);
		if (gam0) gam[jj] = gam0[j];
		else      gam[jj] = (jj == 0) ? 11.0 : (jj <= 2 ? 13.0 : 15.0);
		for (int i = 0; i < qd; i++)
			L[i * p + jj] = L0[i * p0 + j];
		jj++;
	}

	double *mv = (double *)malloc((size_t)qd * sizeof(double));
	for (int i = 0; i < qd; i++)
		mv[i] = m0 ? m0[i] : 1.0;

	/* contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (eq. 2.7) */
	double *r = (double *)malloc((size_t)p * sizeof(double));
	for (int j = 0; j < p; j++)
		r[j] = pow(10.0, -gam[j] / (2.0 * l[j] * N[j]));

	/* restrictive static scaling (eqs. 5.41-5.46), outer vars at |z_k| = r_k */
	double *alpha = (double *)malloc((size_t)p * sizeof(double));
	double *used = (double *)calloc((size_t)qd, sizeof(double));
	int *idx = (int *)malloc((size_t)qd * sizeof(int));
	double *e = (double *)malloc((size_t)qd * sizeof(double));
	for (int j = 0; j < p; j++) {
		int Kj = N[j], lj = l[j];
		int nc = 0;
		for (int i = 0; i < qd; i++) {
			double denom = 1.0 - used[i];
			if (denom <= 0.0) denom = DBL_MIN;
			e[i] = L[i * p + j] / denom;
			if (L[i * p + j] > 0.0) idx[nc++] = i;
		}
		double aj = INFINITY;
		if (nc > 0) {
			/* sort idx by e descending (insertion sort, qd small) */
			for (int a = 1; a < nc; a++) {
				int key = idx[a]; double ek = e[key]; int b = a - 1;
				while (b >= 0 && e[idx[b]] < ek) { idx[b + 1] = idx[b]; b--; }
				idx[b + 1] = key;
			}
			double sumE = 0.0, cummb = 0.0, twolK = 2.0 * lj * Kj;
			for (int n = 0; n < nc; n++) {
				int qi = idx[n];
				sumE += e[qi];
				cummb += mv[qi];
				double cumrho = sumE / (n + 1);
				/* N_{ij} = mbar_n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43) */
				double Nnd = cummb - 1.0;
				for (int k = j + 1; k < p; k++)
					if (L[qi * p + k] != 0.0) Nnd += N[k];
				int Nn = (int)(Nnd + 0.5);
				double an;
				if (Nn <= 0) {
					an = 1.0;
				} else {
					double prod = 1.0;
					for (int ll = 1; ll <= Nn; ll++)
						prod *= (Kj + ll) / (Kj + twolK + ll);
					an = pow(prod, 1.0 / twolK);
				}
				double cand = an / cumrho;
				if (cand < aj) aj = cand;
			}
		}
		if (Z[j] > 0.0) {
			double cand = Kj / Z[j];
			if (cand < aj) aj = cand;
		}
		if (!isfinite(aj)) aj = 1.0;
		alpha[j] = aj;
		for (int i = 0; i < qd; i++)
			used[i] += aj * L[i * p + j] * r[j];
	}

	/* scaled process parameters */
	double *arho0 = (double *)malloc((size_t)p * sizeof(double));
	double *rhoS = (double *)malloc((size_t)qd * p * sizeof(double));
	for (int j = 0; j < p; j++) {
		arho0[j] = alpha[j] * Z[j];
		for (int i = 0; i < qd; i++)
			rhoS[i * p + j] = L[i * p + j] * alpha[j];
	}

	clw_ctx c = { qd, p, N, l, r, arho0, rhoS, mv };
	double complex *w = (double complex *)calloc((size_t)p, sizeof(double complex));
	double gbar = creal(clw_invert(&c, 0, w));

	/* recovery g(N) = prod alpha0_j^{-1} prod alpha_j^{-K_j} gbar (eq. 7.1) */
	double lsum = 0.0;
	for (int j = 0; j < p; j++)
		lsum += arho0[j] - N[j] * log(alpha[j]);
	double lG = log(gbar) + lsum;
	*lGout = lG;
	*Gout = (lG > 709.0) ? INFINITY : exp(lG);

	free(w); free(rhoS); free(arho0); free(e); free(idx); free(used);
	free(alpha); free(r); free(mv); free(L); free(Z); free(gam); free(l);
	free(N); free(keep);
}
