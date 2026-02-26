#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

void freemodel(qnmodel* qn)
{
	if (qn == NULL) return;

	int r, m;

	/* free N array */
	if (qn->N != NULL)
		free(qn->N);

	/* free Z array */
	if (qn->Z != NULL) {
		for (r = 0; r < qn->R; r++)
			mpz_clear(qn->Z[r]);
		free(qn->Z);
	}

	/* free L matrix */
	if (qn->L != NULL) {
		for (m = 0; m < qn->M; m++) {
			if (qn->L[m] != NULL) {
				for (r = 0; r < qn->R; r++)
					mpz_clear(qn->L[m][r]);
				free(qn->L[m]);
			}
		}
		free(qn->L);
	}

	/* free mi array */
	if (qn->mi != NULL)
		free(qn->mi);

	/* free lambda array */
	if (qn->lambda != NULL) {
		for (r = 0; r < qn->R; r++)
			mpq_clear(qn->lambda[r]);
		free(qn->lambda);
	}

	/* free mu matrix */
	if (qn->mu != NULL) {
		for (m = 0; m < qn->M; m++) {
			if (qn->mu[m] != NULL) {
				for (int k = 0; k < qn->Nt; k++)
					mpq_clear(qn->mu[m][k]);
				free(qn->mu[m]);
			}
		}
		free(qn->mu);
	}

	free(qn);
}
