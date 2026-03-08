#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

double lcfsmva_multi(qnmodel* qn, mpq_t *X, mpq_t **Q, mpq_t G, mpq_t **B_out);

int main(int argc, char **argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;

	int r, m;
	qnmodel* qn = NULL;
	/* parse command line arguments */
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	bool exact_output = false;
	bool back_output = false;
	char* model_file = NULL;

	if (argc < 2) {
		printf("USAGE: lcfsmva [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-b|--back] [-d|--exact] [-h|--help] model.qn\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -b, --back       : Print back probabilities (2 rows, R columns)\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--log") == 0) {
			log_output = true;
		} else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ex") == 0) {
			normconst_output = true;
		} else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--nc") == 0) {
			normconst_g_output = true;
		} else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tput") == 0) {
			throughput_output = true;
		} else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qlen") == 0) {
			queue_output = true;
		} else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--back") == 0) {
			back_output = true;
		} else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--exact") == 0) {
			exact_output = true;
		} else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: lcfsmva [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-b|--back] [-d|--exact] [-h|--help] model.qn\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -b, --back       : Print back probabilities (2 rows, R columns)\n");
			printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
			printf("  -h, --help       : Print this help message\n");
			return 0;
		} else if (argv[i][0] == '-') {
			/* Skip unknown options */
			if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0) {
				/* Skip the next argument too (option value) */
				if (i + 1 < argc) i++;
			}
			/* Otherwise just skip this unknown option */
		} else {
			/* Not an option, must be the model file */
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: lcfsmva [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-b|--back] [-d|--exact] [-h|--help] model.qn\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -b, --back       : Print back probabilities (2 rows, R columns)\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}

	qn = (qnmodel*)readmodel(model_file);

	/* Validate model constraints for LCFS MVA */
	if (qn->M != 2) {
		fprintf(stderr, "Error: LCFS MVA requires exactly 2 stations, but model has M=%d\n", qn->M);
		freemodel(qn);
		return -1;
	}
	for (r = 0; r < qn->R; r++) {
		if (mpz_cmp_ui(qn->Z[r], 0) != 0) {
			fprintf(stderr, "Error: LCFS MVA requires Z[%d]=0 (no think times), but Z[%d] is nonzero\n", r + 1, r + 1);
			freemodel(qn);
			return -1;
		}
	}
	for (m = 0; m < qn->M; m++) {
		if (qn->mi[m] != 1) {
			fprintf(stderr, "Error: LCFS MVA requires single-server stations (mi=1), but station %d has mi=%d\n", m + 1, qn->mi[m]);
			freemodel(qn);
			return -1;
		}
	}

	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !back_output && !exact_output) {
		printmodel(qn);
	}

	if (qn->R <= 0 || qn->M <= 0) {
		fprintf(stderr, "Error: Invalid model dimensions R=%d, M=%d\n", qn->R, qn->M);
		freemodel(qn);
		return -1;
	}

	/* Allocate X[R] */
	mpq_t* X = calloc((size_t)qn->R, sizeof(mpq_t));
	if (X == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for X array (%d elements)\n", qn->R);
		freemodel(qn);
		return -1;
	}
	for (r = 0; r < qn->R; r++)
		mpq_init(X[r]);

	/* Allocate Q[2][R] */
	mpq_t** Q = calloc(2, sizeof(mpq_t*));
	if (Q == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for Q array\n");
		for (r = 0; r < qn->R; r++)
			mpq_clear(X[r]);
		free(X);
		freemodel(qn);
		return -1;
	}
	for (m = 0; m < 2; m++) {
		Q[m] = calloc((size_t)qn->R, sizeof(mpq_t));
		if (Q[m] == NULL) {
			fprintf(stderr, "Error: Failed to allocate memory for Q[%d] array\n", m);
			for (r = 0; r < qn->R; r++)
				mpq_clear(X[r]);
			for (int i = 0; i < m; i++) {
				for (r = 0; r < qn->R; r++)
					mpq_clear(Q[i][r]);
				free(Q[i]);
			}
			free(X);
			free(Q);
			freemodel(qn);
			return -1;
		}
		for (r = 0; r < qn->R; r++)
			mpq_init(Q[m][r]);
	}

	/* Allocate B_out[2][R] */
	mpq_t** B_out = calloc(2, sizeof(mpq_t*));
	if (B_out == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for B_out array\n");
		for (r = 0; r < qn->R; r++)
			mpq_clear(X[r]);
		for (m = 0; m < 2; m++) {
			for (r = 0; r < qn->R; r++)
				mpq_clear(Q[m][r]);
			free(Q[m]);
		}
		free(X);
		free(Q);
		freemodel(qn);
		return -1;
	}
	for (m = 0; m < 2; m++) {
		B_out[m] = calloc((size_t)qn->R, sizeof(mpq_t));
		if (B_out[m] == NULL) {
			fprintf(stderr, "Error: Failed to allocate memory for B_out[%d] array\n", m);
			for (r = 0; r < qn->R; r++)
				mpq_clear(X[r]);
			for (int i = 0; i < 2; i++) {
				if (Q[i]) {
					for (r = 0; r < qn->R; r++)
						mpq_clear(Q[i][r]);
					free(Q[i]);
				}
			}
			for (int i = 0; i < m; i++) {
				for (r = 0; r < qn->R; r++)
					mpq_clear(B_out[i][r]);
				free(B_out[i]);
			}
			free(X);
			free(Q);
			free(B_out);
			freemodel(qn);
			return -1;
		}
		for (r = 0; r < qn->R; r++)
			mpq_init(B_out[m][r]);
	}

	mpq_t G;
	mpq_init(G);
	double logG = lcfsmva_multi(qn, X, Q, G, B_out);

	/* Check if solver failed */
	if (logG < 0 && mpq_cmp_ui(G, 0, 1) == 0) {
		for (r = 0; r < qn->R; r++)
			mpq_clear(X[r]);
		for (m = 0; m < 2; m++) {
			for (r = 0; r < qn->R; r++) {
				mpq_clear(Q[m][r]);
				mpq_clear(B_out[m][r]);
			}
			free(Q[m]);
			free(B_out[m]);
		}
		free(Q);
		free(X);
		free(B_out);
		mpq_clear(G);
		freemodel(qn);
		return -1;
	}

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
	} else if (throughput_output) {
		for (r = 1; r <= qn->R; r++) {
			mpf_t xval; mpf_init(xval); mpf_set_q(xval, X[r-1]);
			printf("%.15e\n", mpf_get_d(xval));
			mpf_clear(xval);
		}
	} else if (queue_output) {
		for (m = 1; m <= qn->M; m++) {
			for (r = 1; r <= qn->R; r++) {
				mpf_t qval; mpf_init(qval); mpf_set_q(qval, Q[m-1][r-1]);
				printf("%.15e", mpf_get_d(qval));
				if (r < qn->R) printf(" ");
				mpf_clear(qval);
			}
			printf("\n");
		}
	} else if (back_output) {
		for (m = 0; m < 2; m++) {
			for (r = 0; r < qn->R; r++) {
				mpf_t bval; mpf_init(bval); mpf_set_q(bval, B_out[m][r]);
				printf("%.15e", mpf_get_d(bval));
				if (r < qn->R - 1) printf(" ");
				mpf_clear(bval);
			}
			printf("\n");
		}
	} else if (exact_output) {
		printf("========== Performance Metrics (Exact) ==========\n");
		printf("Normalizing constant:\n");
		gmp_printf("G = %Qd\n", G);
		printf("log(G) = %.15e\n", logG);

		printf("\nX (throughputs):\n");
		for (r = 1; r <= qn->R; r++) {
			gmp_printf("X[%d] = %Qd\n", r, X[r-1]);
		}

		printf("\nQ (mean queue lengths):\n");
		for (m = 1; m <= qn->M; m++) {
			printf("Q[%d] =", m);
			mpq_t total_q; mpq_init(total_q); mpq_set_ui(total_q, 0, 1);
			for (r = 1; r <= qn->R; r++) {
				gmp_printf("\t%Qd", Q[m-1][r-1]);
				mpq_add(total_q, total_q, Q[m-1][r-1]);
			}
			gmp_printf("\t(total: %Qd)\n", total_q);
			mpq_clear(total_q);
		}

		printf("\nB (back probabilities):\n");
		for (m = 0; m < 2; m++) {
			printf("B[%d] =", m + 1);
			for (r = 0; r < qn->R; r++) {
				gmp_printf("\t%Qd", B_out[m][r]);
			}
			printf("\n");
		}
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (LCFS-MVA): %.6f s\n", t1 - t0);
	} else {
		printf("========== Performance Metrics ==========\n");
		printf("Normalizing constant:\n");
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("G = %.15e\n", mpf_get_d(fval));
		printf("log(G) = %.15e\n", logG);
		mpf_clear(fval);
	}

	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !back_output && !exact_output) {
		printf("\nX (throughputs):\n");
		for (r = 1; r <= qn->R; r++) {
			mpf_t xval; mpf_init(xval); mpf_set_q(xval, X[r-1]);
			printf("X[%d] = %.15e\n", r, mpf_get_d(xval));
			mpf_clear(xval);
		}

		printf("\nQ (mean queue lengths):\n");
		for (m = 1; m <= qn->M; m++) {
			printf("Q[%d] =", m);
			mpq_t total_q; mpq_init(total_q); mpq_set_ui(total_q, 0, 1);
			for (r = 1; r <= qn->R; r++) {
				mpf_t qval; mpf_init(qval); mpf_set_q(qval, Q[m-1][r-1]);
				printf("\t%.15e", mpf_get_d(qval));
				mpf_clear(qval);
				mpq_add(total_q, total_q, Q[m-1][r-1]);
			}
			mpf_t tval; mpf_init(tval); mpf_set_q(tval, total_q);
			printf("\t(total: %.15e)\n", mpf_get_d(tval));
			mpf_clear(tval);
			mpq_clear(total_q);
		}

		printf("\nB (back probabilities):\n");
		for (m = 0; m < 2; m++) {
			printf("B[%d] =", m + 1);
			for (r = 0; r < qn->R; r++) {
				mpf_t bval; mpf_init(bval); mpf_set_q(bval, B_out[m][r]);
				printf("\t%.15e", mpf_get_d(bval));
				mpf_clear(bval);
			}
			printf("\n");
		}
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (LCFS-MVA): %.6f s\n", t1 - t0);
	}

	/* Clean up */
	for (r = 0; r < qn->R; r++)
		mpq_clear(X[r]);
	for (m = 0; m < 2; m++) {
		for (r = 0; r < qn->R; r++) {
			mpq_clear(Q[m][r]);
			mpq_clear(B_out[m][r]);
		}
		free(Q[m]);
		free(B_out[m]);
	}
	free(Q);
	free(X);
	free(B_out);
	mpq_clear(G);

	return 0;
}
