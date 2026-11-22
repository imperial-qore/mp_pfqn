#ifndef GLD_H
#define GLD_H

#include <gmp.h>
#include "util.h"

/* Multi-class load-dependent normalizing constant (think time Z included) */
void gld_multi(mpq_t G, mpz_t **L, int *N, mpq_t **mu, mpz_t *Z, int M, int R, int Nt);

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */
void gld_auto_mu(qnmodel *qn);

#endif
