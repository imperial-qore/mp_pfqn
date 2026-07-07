#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "clwld.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

int main(int argc, char **argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;

	bool log_output = false;
	bool normconst_g_output = false;
	char *model_file = NULL;

	if (argc < 2) {
		printf("USAGE: %s [-l] [-g] [-h] model.qn\n", argv[0]);
		printf("  -l : log of normalizing constant\n");
		printf("  -g : normalizing constant as double\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "-g") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-h") == 0) {
			printf("USAGE: %s [-l] [-g] [-h] model.qn\n", argv[0]);
			return 0;
		} else if (argv[i][0] == '-') {
			if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0)
				if (i + 1 < argc) i++;
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: %s [-l] [-g] model.qn\n", argv[0]);
		return -1;
	}

	qnmodel *qn = readmodel(model_file);
	if (!log_output && !normconst_g_output)
		printmodel(qn);

	int qd = qn->M, p = qn->R;
	int ncol = qn->Nt > 0 ? qn->Nt : 1;

	double *L = (double *)malloc((size_t)qd * p * sizeof(double));
	for (int i = 0; i < qd; i++)
		for (int r = 0; r < p; r++)
			L[i * p + r] = mpz_get_d(qn->L[i][r]);
	double *Z = (double *)malloc((size_t)p * sizeof(double));
	for (int r = 0; r < p; r++)
		Z[r] = mpz_get_d(qn->Z[r]);

	/* load-dependent rates S_i(k), k=1..ncol; from MU section or from mi */
	double *mu = (double *)malloc((size_t)qd * ncol * sizeof(double));
	for (int i = 0; i < qd; i++)
		for (int k = 0; k < ncol; k++) {
			if (qn->isLD && qn->mu != NULL)
				mu[i * ncol + k] = mpq_get_d(qn->mu[i][k]);
			else
				mu[i * ncol + k] = (double)MIN(k + 1, qn->mi[i]);
		}

	double G, lG;
	pfqn_clw_lld(qd, p, L, qn->N, Z, mu, ncol, NULL, NULL, &G, &lG);

	if (log_output) {
		printf("%.15e\n", lG);
	} else if (normconst_g_output) {
		printf("%.15e\n", G);
	} else {
		printf("========== CLW-LLD Normalizing Constant ==========\n");
		printf("G = %.15e\n", G);
		printf("log(G) = %.15e\n", lG);
		t1 = CPUTIME;
		printf("Elapsed time (CLW-LLD): %.6f s\n", t1 - t0);
	}

	free(L);
	free(Z);
	free(mu);
	freemodel(qn);
	return 0;
}
