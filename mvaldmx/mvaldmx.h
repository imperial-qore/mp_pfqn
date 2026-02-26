#ifndef MVALDMX_H
#define MVALDMX_H

#include <gmp.h>
#include "util.h"

/*
 * Effective capacity computation for MVA-LD-MX.
 *
 * Inputs:
 *   lambda - arrival rate vector [R], mpq_t
 *   L      - service demand matrix [M][R], mpz_t
 *   mu     - load-dependent rates [M][mu_cols], mpq_t
 *   M      - number of stations
 *   R      - number of classes
 *   Nt     - total closed population
 *   mu_cols - number of columns in mu (>= Nt+1)
 *
 * Outputs (caller must allocate/init):
 *   EC     - effective capacity [M][Nt], indexed EC[m][n-1] for n=1..Nt
 *   E      - E values [M][Nt+1], indexed E[m][n] for n=0..Nt
 *   Eprime - E' values [M][Nt+1], indexed Eprime[m][n] for n=0..Nt
 */
void mvaldmx_ec(mpq_t **EC, mpq_t **E, mpq_t **Eprime,
                mpq_t *lambda, mpz_t **L, mpq_t **mu,
                int M, int R, int Nt, int mu_cols);

/*
 * MVA-LD-MX solver: load-dependent MVA for mixed open/closed networks.
 *
 * Inputs:
 *   qn - queueing network model (with hasOpen, lambda, isLD, mu fields)
 *
 * Outputs (caller must allocate/init):
 *   X  - throughputs [R], mpq_t
 *   Q  - queue lengths [M][R], mpq_t
 *
 * Returns: log of normalizing constant (as double), or -1 on error.
 */
double mvaldmx_solve(qnmodel *qn, mpq_t *X, mpq_t **Q);

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) for k=0..Nt-1 */
void mvaldmx_auto_mu(qnmodel *qn);

#endif
