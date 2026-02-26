#ifndef GLD_H
#define GLD_H

#include <gmp.h>
#include "util.h"

/* Single-class load-dependent DP normalizing constant */
void gld_single(mpq_t G, mpz_t *L, int N, mpq_t **mu, int M);

/* Multi-class load-dependent recursion */
void gld_multi(mpq_t G, mpz_t **L, int *N, mpq_t **mu, int M, int R, int Nt);

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */
void gld_auto_mu(qnmodel *qn);

#endif
