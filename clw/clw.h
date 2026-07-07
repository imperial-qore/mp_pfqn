#ifndef CLW_H
#define CLW_H

#include "util.h"

/*
 * Choudhury-Leung-Whitt normalizing constant by numerical inversion of the
 * generating function (JACM 42(5):935-970, 1995).
 *
 * Computes g(N) of a multichain closed product-form network with single-server
 * and (optionally) infinite-server queues by numerically inverting its
 * p-dimensional generating function (eq. 4.5)
 *
 *   G(z) = exp(sum_j rho_{j0} z_j) / prod_i (1 - sum_j rho_{ji} z_j)^{m_i}.
 *
 * All arrays are row-major doubles. This is a floating-point method: G and lG
 * are returned as doubles (about 7-9 significant digits), not exact rationals.
 *
 * @param qd   number of distinct single-server queues (rows of L).
 * @param p    number of closed chains (columns of L).
 * @param L    (qd*p) relative traffic intensities, L[i*p+j] = rho_{ji}.
 * @param N    (p) population vector K.
 * @param Z    (p) aggregate infinite-server intensities rho_{j0}, or NULL for 0.
 * @param m    (qd) queue multiplicities m_i, or NULL for all ones.
 * @param l    (p) inner lattice parameters l_j, or NULL for the CLW defaults.
 * @param gam  (p) aliasing parameters gamma_j, or NULL for the CLW defaults.
 * @param G    output normalizing constant g(N) (Inf on overflow).
 * @param lG   output natural logarithm of g(N) (always finite).
 */
void pfqn_clw(int qd, int p, const double *L, const int *N, const double *Z,
              const double *m, const int *l, const double *gam,
              double *G, double *lG);

#endif
