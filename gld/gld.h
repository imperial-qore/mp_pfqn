#ifndef GLD_H
#define GLD_H

#include <gmp.h>
#include "util.h"

/* Multi-class load-dependent normalizing constant (think time Z included) */
void gld_multi(mpq_t G, mpz_t **L, int *N, mpq_t **mu, mpz_t *Z, int M, int R, int Nt);

/* P = prod_r (N_r+1), the size of the population box */
long int gld_boxsize(int *N, int R);

/* G(n) for all n <= N into Gv[P] (mpq_init'ed by the caller), mixed-radix
 * index sum_r n_r prod_{s<r}(N_s+1); station skip (-1: none) left out */
void gld_table(mpq_t *Gv, mpz_t **L, int *N, mpq_t **mu, mpz_t *Z, int M, int R, int Nt, int skip);

/* balance factor Y(k), k <= N, of a station (mu_m == NULL: the delay, dem = Z) */
void gld_station_factor(mpq_t *Y, mpz_t *dem, mpq_t *mu_m, int *N, int R, int Nt);

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */

#endif
