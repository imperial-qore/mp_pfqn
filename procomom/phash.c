#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "procomom.h"

/* phash - map (Dn combination, station index) to basis index (0-based)
   comb: an R-length Dn combination vector
   i: station index, 1..M
   Returns: 0-based index into pk vector of size cardG * M
            or -1 if combination not found */
int phash(combsrep* Dn, int* comb, int i)
{
	long int pos = int_matmatchrow(Dn->combs, Dn->card, Dn->n, comb);
	if (pos < 0) {
		return -1;
	}
	return (int)(pos * Dn->k + (i - 1));
}
