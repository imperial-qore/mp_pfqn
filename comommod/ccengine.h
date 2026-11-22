#ifndef CCENGINE_H
#define CCENGINE_H

/* API of the multi-modular CoMoM engine (comommod/engine.c), the CoMoM
 * counterpart of mommod/mmengine.h.  bin/comom links it for tier 4 of its
 * singular-model fallback. */

#include <stddef.h>
#include <gmp.h>
#include "util.h"

/* Prime count per batch (OpenMP threads), prime size in bits, witness
 * primes and bounded reconstruction (-B); denbits > 0 fixes the
 * denominator bound at 2^denbits, 0 takes a third of the modulus. */
void comommod_config(int nthreads, int primebits, int nwitness, int use_bounded, int denbits);
/* Primes used, unlucky primes discarded, and modulus bits of the last run. */
void comommod_stats(long* primes, long* unlucky, unsigned long* modbits);
/* Exact exported constants of q (layout of util/singular.h) through every
 * tier: 0 solved, 1 every tier failed, 2 hard error already reported. */
int comommod_solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz);
/* Tier 4 alone: interpolation of the exact limit of a perturbation, per
 * prime.  Same contract; usable as sg_ops.interp (ctx is ignored). */
int comommod_interp(qnmodel* q, mpq_t* out, void* ctx);

#endif
