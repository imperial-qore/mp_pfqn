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
#include "replicas.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

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
		printf("USAGE: %s [-l] [-e] [-g] [-P] [-h] model.qn\n", argv[0]);
		printf("  -l : log of normalizing constant\n");
		printf("  -e : exact normalizing constant (num/den)\n");
		printf("  -g : normalizing constant as double\n");
		printf("  -P : state probabilities at the queueing station (was -q)\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0)
			normconst_output = true;
		else if (strcmp(argv[i], "-g") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-P") == 0)
			prob_output = true;
		else if (strcmp(argv[i], "-q") == 0) {
			/* -q means mean queue lengths in every other solver */
			fprintf(stderr, "comomld: the state probabilities moved from -q to -P; mean queue lengths (-q) are not available\n");
			return 1;
		}
		else if (strcmp(argv[i], "-h") == 0) {
			printf("USAGE: %s [-l] [-e] [-g] [-P] [-h] model.qn\n", argv[0]);
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
		printf("USAGE: %s [-l] [-e] [-g] [-P] model.qn\n", argv[0]);
		return -1;
	}

	qnmodel *qn = readmodel(model_file);
	qn_require_closed(qn, "comomld", 1);
	if (!log_output && !normconst_output && !normconst_g_output && !prob_output)
		printmodel(qn);

	/* mi copies of a station as one load-dependent station (util/replicas.h) */
	if (qn_replica_mu(qn))
		return 1;

	int Nt = qn->Nt;

	mpq_t G;
	mpq_init(G);

	/* Allocate probability vector */
	mpq_t *prob = (mpq_t *)malloc((Nt + 1) * sizeof(mpq_t));
	for (int k = 0; k <= Nt; k++)
		mpq_init(prob[k]);

	int status = comomrm_ld(G, prob, qn);

	if (status != 0) {
		/* Model unsupported by CoMoM-LD (comomrm_ld has said why): G is
		 * undefined.  Fail like every other solver does outside its domain,
		 * rather than printing NaN with exit status 0 (which in turn had
		 * replaced a silent G = 0). */
		fprintf(stderr, "comomld: model outside the CoMoM-LD domain; no result\n");
		for (int k = 0; k <= Nt; k++)
			mpq_clear(prob[k]);
		free(prob);
		mpq_clear(G);
		freemodel(qn);
		return 1;
	}

	/* Compute logG */
	/* mpq_log: log(mpf_get_d(G)) is inf once G exceeds ~1e308 */
	double logG = mpq_sgn(G) > 0 ? mpq_log(G) : log(0.0);

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
