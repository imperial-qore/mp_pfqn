#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <zpla.h>
#include "mommod.h"

/* Collect from one modular image the residues of every quantity the
 * solver reports, in the layout
 *
 *   res[0]                        = G
 *   res[1 + (r-1)]                = G_r
 *   res[1 + R + (k-1)*R + (r-1)]  = G^{+k}_r
 *
 * The raw constants are exported rather than the performance measures
 * X_r = G_r/G and Q_kr = m_k rho_kr G^{+k}_r/G.  Forming the ratio inside
 * the ring would be cheaper per prime, but the ratio is a rational whose
 * numerator and denominator are each about as long as G, so rational
 * reconstruction of it needs roughly four times the modulus that the
 * integers themselves need, hence four times as many primes.  The
 * division is done once, in exact arithmetic, after reconstruction. */
void exportres(qnmodel* qnm, zp_vec_t G, zp_vec_t Gk, zp_t* res)
{
	int r,k;
	int R = qnm->R;

	res[0] = G[0];
	for (r=1;r<=R;r++) res[r] = G[r];
	for (k=1;k<=qnm->M;k++)
		for (r=1;r<=R;r++)
			res[1 + R + (k-1)*R + (r-1)] = Gk[(k-1)*(R+1)+r];
}
