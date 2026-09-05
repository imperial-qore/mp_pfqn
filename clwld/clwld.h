#ifndef CLWLD_H
#define CLWLD_H

#include "util.h"

/*
 * Choudhury-Leung-Whitt normalizing constant by numerical inversion of the
 * generating function, extended to limited load-dependent (LLD) stations via
 * the per-center transforms of Bertozzi and McKenna (SIAM Review 35(2):239-268,
 * 1993).
 *
 * The generating function is (Bertozzi-McKenna eqs. 2.17/2.23)
 *
 *   G(z) = exp(sum_j rho_{j0} z_j) prod_i F_i(sum_j rho_{ji} z_j),
 *
 * with F_i(x) = sum_{n>=0} x^n / prod_{k=1}^n S_i(k) the transform of the
 * station factor of queue i and S_i(k) = mu(i,k) its load-dependent rate. For
 * an LLD queue (S_i(k) = c_i for k >= l_i) F_i is rational with a simple pole
 * at x = c_i. Multiserver (S_i(k) = min(k,c_i)) and load-independent
 * (F_i = 1/(1-x)) queues are special cases.
 *
 * Floating-point method: G and lG are returned as doubles, not exact rationals.
 *
 * @param qd   number of queues (rows of L and mu).
 * @param p    number of closed chains (columns of L).
 * @param L    (qd*p) relative traffic intensities, L[i*p+j] = rho_{ji}.
 * @param N    (p) population vector K.
 * @param Z    (p) aggregate infinite-server intensities rho_{j0}, or NULL for 0.
 * @param mu   (qd*ncol) load-dependent rates, mu[i*ncol+(k-1)] = S_i(k).
 * @param ncol number of load-dependent columns (>= sum(N)); the last column
 *             is the LLD constant rate c_i.
 * @param l    (p) inner lattice parameters l_j, or NULL for the CLW defaults.
 * @param gam  (p) aliasing parameters gamma_j, or NULL for the CLW defaults.
 * @param G    output normalizing constant g(N) (Inf on overflow).
 * @param lG   output natural logarithm of g(N) (always finite).
 */
void pfqn_clw_lld(int qd, int p, const double *L, const int *N, const double *Z,
                  const double *mu, int ncol, const int *l, const double *gam,
                  double *G, double *lG);

#endif
