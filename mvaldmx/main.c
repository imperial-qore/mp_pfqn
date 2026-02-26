#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "mvaldmx.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

int main(int argc, char **argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;

	int r, m;
	qnmodel *qn = NULL;

	/* Parse command line arguments */
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	bool exact_output = false;
	char *model_file = NULL;

	if (argc < 2) {
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--log") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ex") == 0)
			normconst_output = true;
		else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--nc") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tput") == 0)
			throughput_output = true;
		else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qlen") == 0)
			queue_output = true;
		else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--exact") == 0)
			exact_output = true;
		else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
			return 0;
		} else if (argv[i][0] == '-') {
			/* Skip unknown options */
			if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0)
				if (i + 1 < argc) i++;
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] model.qn\n", argv[0]);
		return -1;
	}

	qn = readmodel(model_file);
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !exact_output)
		printmodel(qn);

	if (qn->R <= 0 || qn->M <= 0) {
		fprintf(stderr, "Error: Invalid model dimensions R=%d, M=%d\n", qn->R, qn->M);
		return -1;
	}

	if (qn->Nt <= 0) {
		fprintf(stderr, "Error: No closed classes found (Nt=%d)\n", qn->Nt);
		return -1;
	}

	/* Auto-generate mu from mi if no MU section */
	if (!qn->isLD)
		mvaldmx_auto_mu(qn);

	/* Allocate output arrays */
	mpq_t *X = (mpq_t *)calloc(qn->R, sizeof(mpq_t));
	mpq_t **Q = (mpq_t **)calloc(qn->M, sizeof(mpq_t *));
	for (r = 0; r < qn->R; r++)
		mpq_init(X[r]);
	for (m = 0; m < qn->M; m++) {
		Q[m] = (mpq_t *)calloc(qn->R, sizeof(mpq_t));
		for (r = 0; r < qn->R; r++)
			mpq_init(Q[m][r]);
	}

	/* Run solver */
	double logG = mvaldmx_solve(qn, X, Q);

	/* Output results */
	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		/* MVA-LD-MX does not compute a normalizing constant directly */
		printf("0\n1\n");
	} else if (normconst_g_output) {
		printf("%.15e\n", 0.0);
	} else if (throughput_output) {
		for (r = 0; r < qn->R; r++) {
			mpf_t xval;
			mpf_init(xval);
			mpf_set_q(xval, X[r]);
			printf("%.15e\n", mpf_get_d(xval));
			mpf_clear(xval);
		}
	} else if (queue_output) {
		for (m = 0; m < qn->M; m++) {
			for (r = 0; r < qn->R; r++) {
				mpf_t qval;
				mpf_init(qval);
				mpf_set_q(qval, Q[m][r]);
				printf("%.15e", mpf_get_d(qval));
				if (r < qn->R - 1) printf(" ");
				mpf_clear(qval);
			}
			printf("\n");
		}
	} else if (exact_output) {
		printf("========== Performance Metrics (Exact) ==========\n");

		printf("\nX (throughputs):\n");
		for (r = 0; r < qn->R; r++)
			gmp_printf("X[%d] = %Qd\n", r + 1, X[r]);

		printf("\nQ (mean queue lengths):\n");
		for (m = 0; m < qn->M; m++) {
			printf("Q[%d] =", m + 1);
			mpq_t total_q;
			mpq_init(total_q);
			for (r = 0; r < qn->R; r++) {
				gmp_printf("\t%Qd", Q[m][r]);
				mpq_add(total_q, total_q, Q[m][r]);
			}
			gmp_printf("\t(total: %Qd)\n", total_q);
			mpq_clear(total_q);
		}
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (MVALDMX): %.6f s\n", t1 - t0);
	} else {
		printf("========== Performance Metrics ==========\n");

		printf("\nX (throughputs):\n");
		for (r = 0; r < qn->R; r++) {
			mpf_t xval;
			mpf_init(xval);
			mpf_set_q(xval, X[r]);
			printf("X[%d] = %.15e\n", r + 1, mpf_get_d(xval));
			mpf_clear(xval);
		}

		printf("\nQ (mean queue lengths):\n");
		for (m = 0; m < qn->M; m++) {
			printf("Q[%d] =", m + 1);
			mpq_t total_q;
			mpq_init(total_q);
			for (r = 0; r < qn->R; r++) {
				mpf_t qval;
				mpf_init(qval);
				mpf_set_q(qval, Q[m][r]);
				printf("\t%.15e", mpf_get_d(qval));
				mpf_clear(qval);
				mpq_add(total_q, total_q, Q[m][r]);
			}
			mpf_t tval;
			mpf_init(tval);
			mpf_set_q(tval, total_q);
			printf("\t(total: %.15e)\n", mpf_get_d(tval));
			mpf_clear(tval);
			mpq_clear(total_q);
		}
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (MVALDMX): %.6f s\n", t1 - t0);
	}

	/* Cleanup */
	for (r = 0; r < qn->R; r++)
		mpq_clear(X[r]);
	for (m = 0; m < qn->M; m++) {
		for (r = 0; r < qn->R; r++)
			mpq_clear(Q[m][r]);
		free(Q[m]);
	}
	free(Q);
	free(X);

	freemodel(qn);
	return 0;
}
