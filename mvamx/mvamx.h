#ifndef MVAMX_H
#define MVAMX_H

#include <gmp.h>
#include "util.h"

/**
 * MVA-MX: Mean Value Analysis for mixed open/closed queueing networks.
 *
 * Computes throughputs, queue lengths, and the log normalizing constant
 * for a product-form network that may contain both open classes (N[r]=-1,
 * arrival rate lambda[r]) and closed classes (N[r]>=0).
 *
 * @param qn      Model (must have hasOpen=1 and lambda[] set for open classes)
 * @param X       Output: throughputs [R] (mpq_t, pre-initialized)
 * @param Q       Output: queue lengths [M][R] (mpq_t, pre-initialized)
 * @param G       Output: normalizing constant of the closed sub-model (mpq_t, pre-initialized)
 * @return        log(G) of the closed sub-model, or 0.0 if no closed classes
 */
double mvamx(qnmodel* qn, mpq_t *X, mpq_t **Q, mpq_t G);

#endif /* MVAMX_H */
