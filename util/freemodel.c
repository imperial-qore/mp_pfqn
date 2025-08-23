#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

void freemodel(qnmodel* qn)
{
	if (qn == NULL) return;
	
	int r, m;
	
	// Free N array
	if (qn->N != NULL) {
		for (r = 0; r < qn->R; r++) {
			mpz_clear(qn->N[r]);
		}
		free(qn->N);
	}
	
	// Free Z array
	if (qn->Z != NULL) {
		for (r = 0; r < qn->R; r++) {
			mpz_clear(qn->Z[r]);
		}
		free(qn->Z);
	}
	
	// Free L matrix
	if (qn->L != NULL) {
		for (m = 0; m < qn->M; m++) {
			if (qn->L[m] != NULL) {
				for (r = 0; r < qn->R; r++) {
					mpz_clear(qn->L[m][r]);
				}
				free(qn->L[m]);
			}
		}
		free(qn->L);
	}
	
	// Free mi array
	if (qn->mi != NULL) {
		for (m = 0; m < qn->M; m++) {
			mpz_clear(qn->mi[m]);
		}
		free(qn->mi);
	}
	
	// Free the structure itself
	free(qn);
}