#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "procomom.h"

/* pexact - compute initial pk for empty network (Ncur = [0,...,0], n = 0)
   For the empty network, only the zero combination is valid (all others
   have N-Dn[d] with negative entries). At n=0, P(0 customers) = 1 for
   all stations.
   Sets pk[phash(zero_comb, k)] = 1 for every component, rest = 0. */
void pexact(mpq_vec_t pk, combsrep* Dn, int M)
{
	(void) M;
	int basisSize = Dn->card * Dn->stride;
	int i, k;

	/* Zero everything */
	for (i = 0; i < basisSize; i++)
		mpq_set_si(pk[i], 0, 1);

	/* The zero combination is the first row in Dn after sortbynnzpos
	   (all zeros sort first). Set pk = 1 for all stations at this
	   combination. */
	int* zero_comb = (int*)int_vec(Dn->n, 0);
	for (k = 1; k <= Dn->stride; k++) {
		int idx = phash(Dn, zero_comb, k);
		if (idx >= 0 && idx < basisSize)
			mpq_set_si(pk[idx], 1, 1);
	}
	free(zero_comb);
}
