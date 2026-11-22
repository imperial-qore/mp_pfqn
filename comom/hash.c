#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"

int hash(combsrep* Dn, int *comb, int i)
{
	// Following JCoMoM's systematic indexing: population_position * (M + 1) + multiplicity
	// where i represents the multiplicity (0 for G, 1..M for G^k)
	// i=0 means no queue added (for G values)
	// i=1..M means queue i added (for G^k values)
	
	int population_position = int_matmatchrow(Dn->combs, Dn->card, Dn->n, comb);
	if (population_position < 0) {
		fprintf(stderr, "Error: hash() - combination not found in basis\n");
		return -1;
	}
	
	// The basis is organized as:
	// - First cardGk = cardG * M entries are G^k values
	// - Next cardG entries are G values
	// So for each population vector, we have M entries for G^k followed by 1 entry for G
	
	if (i == 0) {
		// G value (no queue added) - stored after all G^k values
		return Dn->card * Dn->k + population_position + 1;
	} else {
		// G^k value (queue i added) - stored in first part
		// position = population_index * M + (queue_index - 1)
		return population_position * Dn->k + (i - 1) + 1;
	}
}


