#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "gld.h"
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
	char *model_file = NULL;

	if (argc < 2) {
		printf("USAGE: %s [-l] [-e] [-g] [-h] model.qn\n", argv[0]);
		printf("  -l : log of normalizing constant\n");
		printf("  -e : exact normalizing constant (num/den)\n");
		printf("  -g : normalizing constant as double\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0)
			normconst_output = true;
		else if (strcmp(argv[i], "-g") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-h") == 0) {
			printf("USAGE: %s [-l] [-e] [-g] [-h] model.qn\n", argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0) {
			/* perturbation options of mom/comom: accepted and ignored */
			if (i + 1 < argc) i++;
		} else if (argv[i][0] == '-') {
			fprintf(stderr, "%s: unsupported option '%s' (gld computes the normalizing constant only: use -l, -e or -g)\n", argv[0], argv[i]);
			return 1;
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: %s [-l] [-e] [-g] model.qn\n", argv[0]);
		return -1;
	}

	qnmodel *qn = readmodel(model_file);
	qn_require_closed(qn, "gld", 1);
	if (!log_output && !normconst_output && !normconst_g_output)
		printmodel(qn);

	/* mi copies of a station as one load-dependent station (util/replicas.h) */
	if (qn_replica_mu(qn))
		return 1;

	mpq_t G;
	mpq_init(G);

	gld_multi(G, qn->L, qn->N, qn->mu, qn->Z, qn->M, qn->R, qn->Nt);

	/* Compute logG (mpq_log does not overflow when G > DBL_MAX) */
	double logG = mpq_log(G);

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
	} else {
		printf("========== GLD Normalizing Constant ==========\n");
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("G = %.15e\n", mpf_get_d(fval));
		printf("log(G) = %.15e\n", logG);
		mpf_clear(fval);
		t1 = CPUTIME;
		printf("Elapsed time (GLD): %.6f s\n", t1 - t0);
	}

	mpq_clear(G);
	freemodel(qn);
	return 0;
}
