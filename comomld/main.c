#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "comomld.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

/* Auto-generate mu from mi: mu[m][k] = min(k+1, mi[m]) */
static void comomld_auto_mu(qnmodel *qn)
{
	int m, k;
	if (qn->mu != NULL || qn->Nt <= 0)
		return;
	qn->mu = (mpq_t **)malloc(qn->M * sizeof(mpq_t *));
	for (m = 0; m < qn->M; m++) {
		qn->mu[m] = (mpq_t *)malloc(qn->Nt * sizeof(mpq_t));
		for (k = 0; k < qn->Nt; k++) {
			mpq_init(qn->mu[m][k]);
			mpq_set_ui(qn->mu[m][k], MIN(k + 1, qn->mi[m]), 1);
		}
	}
	qn->isLD = 1;
}

int main(int argc, char **argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;

	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool prob_output = false;
	char *model_file = NULL;

	if (argc < 2) {
		printf("USAGE: %s [-l] [-e] [-g] [-q] [-h] model.qn\n", argv[0]);
		printf("  -l : log of normalizing constant\n");
		printf("  -e : exact normalizing constant (num/den)\n");
		printf("  -g : normalizing constant as double\n");
		printf("  -q : state probabilities at queueing station\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0)
			normconst_output = true;
		else if (strcmp(argv[i], "-g") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-q") == 0)
			prob_output = true;
		else if (strcmp(argv[i], "-h") == 0) {
			printf("USAGE: %s [-l] [-e] [-g] [-q] [-h] model.qn\n", argv[0]);
			return 0;
		} else if (argv[i][0] == '-') {
			/* skip unknown options with arguments */
			if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0)
				if (i + 1 < argc) i++;
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: %s [-l] [-e] [-g] [-q] model.qn\n", argv[0]);
		return -1;
	}

	qnmodel *qn = readmodel(model_file);
	if (!log_output && !normconst_output && !normconst_g_output && !prob_output)
		printmodel(qn);

	/* Auto-generate mu from mi if no MU section */
	if (!qn->isLD)
		comomld_auto_mu(qn);

	int Nt = qn->Nt;

	mpq_t G;
	mpq_init(G);

	/* Allocate probability vector */
	mpq_t *prob = (mpq_t *)malloc((Nt + 1) * sizeof(mpq_t));
	for (int k = 0; k <= Nt; k++)
		mpq_init(prob[k]);

	comomrm_ld(G, prob, qn);

	/* Compute logG */
	mpf_t G_mpf;
	mpf_init(G_mpf);
	mpf_set_q(G_mpf, G);
	double logG = log(mpf_get_d(G_mpf));
	mpf_clear(G_mpf);

	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		mpq_canonicalize(G);
		mpz_t num, den;
		mpz_init(num);
		mpz_init(den);
		mpq_get_num(num, G);
		mpq_get_den(den, G);
		gmp_printf("%Zd\n", num);
		gmp_printf("%Zd\n", den);
		mpz_clear(num);
		mpz_clear(den);
	} else if (normconst_g_output) {
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("%.15e\n", mpf_get_d(fval));
		mpf_clear(fval);
	} else if (prob_output) {
		for (int k = 0; k <= Nt; k++) {
			mpq_canonicalize(prob[k]);
			mpf_t fval;
			mpf_init(fval);
			mpf_set_q(fval, prob[k]);
			printf("%.15e\n", mpf_get_d(fval));
			mpf_clear(fval);
		}
	} else {
		printf("========== CoMoM-LD Normalizing Constant ==========\n");
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("G = %.15e\n", mpf_get_d(fval));
		printf("log(G) = %.15e\n", logG);
		mpf_clear(fval);

		printf("\n========== State Probabilities ==========\n");
		for (int k = 0; k <= Nt; k++) {
			mpq_canonicalize(prob[k]);
			mpf_t pval;
			mpf_init(pval);
			mpf_set_q(pval, prob[k]);
			printf("P(%d) = %.15e\n", k, mpf_get_d(pval));
			mpf_clear(pval);
		}

		t1 = CPUTIME;
		printf("\nElapsed time (CoMoM-LD): %.6f s\n", t1 - t0);
	}

	/* Cleanup */
	for (int k = 0; k <= Nt; k++)
		mpq_clear(prob[k]);
	free(prob);
	mpq_clear(G);
	freemodel(qn);
	return 0;
}
