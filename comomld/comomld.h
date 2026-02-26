#ifndef COMOMLD_H
#define COMOMLD_H

#include <gmp.h>
#include "util.h"

/*
 * CoMoM-LD: Class-Oriented Method of Moments for repairman models
 * with load-dependent service rates.
 *
 * Requires exactly 1 queueing station after infinite-server absorption.
 * Returns normalizing constant G and state probabilities prob[0..Nt].
 *
 * prob[k] = P(k jobs at the queueing station), k=0..Nt
 */
void comomrm_ld(mpq_t G, mpq_t *prob, qnmodel *qn);

#endif
