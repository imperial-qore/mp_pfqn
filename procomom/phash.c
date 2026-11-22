#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "procomom.h"

/* phash - map (Dn combination, component index) to basis index (0-based)
   comb: an R-length Dn combination vector
   i: component index, 1..M+1.  i=1 is the base component G(N-Dn);
      i=k+1 is the component with one extra copy of station k, k=1..M.
   Returns: 0-based index into pk vector of size cardG * (M+1)
            or -1 if combination not found

   The block stride is Dn->stride.  comom's own basis carries the base plus
   ALL M stations (comom/hash.c); this file carries the base plus stations
   1..M-1, dropping station M, because with a single-server reference station
   that component has weight mi_M - 1 = 0 everywhere -- including it would
   leave an unconstrained column and a singular A^T A.  When the reference
   station IS replicated the term is nonzero, and main.c widens the stride to
   M+1 so it has somewhere to live. */
int phash(combsrep* Dn, int* comb, int i)
{
	long int pos = int_matmatchrow(Dn->combs, Dn->card, Dn->n, comb);
	if (pos < 0) {
		return -1;
	}
	return (int)(pos * Dn->stride + (i - 1));
}
