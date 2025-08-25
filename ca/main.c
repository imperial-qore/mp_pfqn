#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"
#include "popcycle.h"

/* Define global variables from profiling.h */
int AORSCTR;
int MULCTR;
int DIVCTR;
double AORSTIME;
double MULTIME;
double DIVTIME;
double t0,t1;
struct rusage ruse;

/* Define global variables from mom.h */
qnmodel* qnm;
bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;

void convolution_multi_exact(mpq_t *g, mpf_t *X, mpf_t **Q);

void printcompact(int*n,int R)
{
	int s;
	fprintf(stdout,"n=(");
	for (s=0;s<R-1;s++)
	    fprintf(stdout,"%d,",n[s]);
	fprintf(stdout,"%d) \n",n[R-1]);
	fflush(stdout);
}


int main(int argc, char**argv)
{
	t0=CPUTIME;

	/* parse command line arguments */
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	bool exact_output = false;
	char* model_file = NULL;
	
	if(argc < 2)
	{
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log      : Print only log of normalizing constant as double\n");
		printf("  -e, --ex       : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact    : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}
	
	for(int i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--log") == 0) {
			log_output = true;
		} else if(strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ex") == 0) {
			normconst_output = true;
		} else if(strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--nc") == 0) {
			normconst_g_output = true;
		} else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tput") == 0) {
			throughput_output = true;
		} else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qlen") == 0) {
			queue_output = true;
		} else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--exact") == 0) {
			exact_output = true;
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log      : Print only log of normalizing constant as double\n");
			printf("  -e, --ex       : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --exact    : Print all performance metrics in full exact precision (integer or rational)\n");
			printf("  -h, --help       : Print this help message\n");
			return 0;
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-l|--log] [-n|--norm-const] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -n, --norm-const : Print exact normalizing constant numerator and denominator\n");
		return -1;
	}
	
	qnm=(qnmodel*)readmodel(model_file);
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !exact_output) {
		printmodel(qnm);
	}
	
	if (qnm->R <= 0 || qnm->M <= 0) {
		fprintf(stderr, "Error: Invalid model dimensions R=%d, M=%d\n", qnm->R, qnm->M);
		return -1;
	}
	
	mpq_t* G=(mpq_t*) calloc((size_t)(qnm->M+1),sizeof(mpq_t));
	if (G == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for G array (%d elements)\n", qnm->M+1);
		return -1;
	}
	
	mpf_t* X=(mpf_t*) calloc((size_t)qnm->R,sizeof(mpf_t));
	if (X == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for X array (%d elements)\n", qnm->R);
		free(G);
		return -1;
	}
	
	mpf_t** Q=(mpf_t**) calloc((size_t)qnm->M,sizeof(mpf_t*));
	if (Q == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for Q array (%d elements)\n", qnm->M);
		free(G);
		free(X);
		return -1;
	}
	for (int i=0;i<=qnm->M;i++) mpq_init(G[i]);
	for (int i=0;i<qnm->R;i++) mpf_init(X[i]);
	for (int i=0;i<qnm->M;i++) {
		Q[i] = (mpf_t*) calloc((size_t)qnm->R,sizeof(mpf_t));
		if (Q[i] == NULL) {
			fprintf(stderr, "Error: Failed to allocate memory for Q[%d] array (%d elements)\n", i, qnm->R);
			// Clean up previously allocated memory
			for (int j=0;j<=qnm->M;j++) mpq_clear(G[j]);
			for (int j=0;j<qnm->R;j++) mpf_clear(X[j]);
			for (int j=0;j<i;j++) {
				for (int r=0;r<qnm->R;r++) mpf_clear(Q[j][r]);
				free(Q[j]);
			}
			free(G);
			free(X);
			free(Q);
			return -1;
		}
		for (int r=0;r<qnm->R;r++) mpf_init(Q[i][r]);
	}
	
	convolution_multi_exact(G,X,Q);
	t1=CPUTIME;
	
	// Check if convolution failed (G[0] would be uninitialized/zero)
	if (mpq_cmp_ui(G[0], 0, 1) == 0) {
		// Clean up and exit gracefully
		for (int i=0;i<=qnm->M;i++) mpq_clear(G[i]);
		for (int i=0;i<qnm->R;i++) mpf_clear(X[i]);
		for (int i=0;i<qnm->M;i++) {
			for (int r=0;r<qnm->R;r++) mpf_clear(Q[i][r]);
			free(Q[i]);
		}
		free(G);
		free(X);
		free(Q);
		return -1;
	}
	
	mpf_t fval;
	mpf_init(fval);
	mpf_set_q(fval, G[0]);
	double logG = log(mpf_get_d(fval));
	
	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		// Print exact numerator and denominator
		mpz_t num, den;
		mpz_init(num);
		mpz_init(den);
		mpq_get_num(num, G[0]);
		mpq_get_den(den, G[0]);
		gmp_printf("%Zd\n", num);
		gmp_printf("%Zd\n", den);
		mpz_clear(num);
		mpz_clear(den);
	} else if (normconst_g_output) {
		// Print normalizing constant as double
		printf("%.15e\n", mpf_get_d(fval));
	} else if (throughput_output) {
		// Print only throughputs, one per row
		for (int r = 0; r < qnm->R; r++) {
			printf("%.15e\n", mpf_get_d(X[r]));
		}
	} else if (queue_output) {
		// Print queue lengths, all classes for same queue on same row
		for (int i = 0; i < qnm->M; i++) {
			for (int r = 0; r < qnm->R; r++) {
				printf("%.15e", mpf_get_d(Q[i][r]));
				if (r < qnm->R - 1) printf(" ");
			}
			printf("\n");
		}
	} else if (exact_output) {
		// Print all performance metrics in exact rational form
		printf("========== Performance Metrics (Exact) ==========\n");
		gmp_printf("G = %Qd\n", G[0]);
		printf("log(G) = %.15e\n", logG);
		
		printf("\nX (throughputs):\n");
		mpq_t x_exact;
		mpq_init(x_exact);
		for (int r = 0; r < qnm->R; r++) {
			mpq_set_f(x_exact, X[r]);
			gmp_printf("X[%d] = %Qd\n", r, x_exact);
		}
		mpq_clear(x_exact);
		
		printf("\nQ (mean queue lengths):\n");
		mpq_t q_exact, q_total;
		mpq_init(q_exact);
		mpq_init(q_total);
		for (int i = 0; i < qnm->M; i++) {
			printf("Q[%d] =", i+1);
			mpq_set_ui(q_total, 0, 1);
			for (int r = 0; r < qnm->R; r++) {
				mpq_set_f(q_exact, Q[i][r]);
				gmp_printf("\t%Qd", q_exact);
				mpq_add(q_total, q_total, q_exact);
			}
			gmp_printf("\t(total: %Qd)\n", q_total);
		}
		mpq_clear(q_exact);
		mpq_clear(q_total);
		
		printf("=========================================\n");
		printf("Elapsed time (CA): %g s\n", t1-t0);
	} else {
		printf("========== Performance Metrics ==========\n");
		printf("G = %.15e\n", mpf_get_d(fval));
		printf("log(G) = %.15e\n", logG);
		
		printf("\nX (throughputs):\n");
		for (int r = 0; r < qnm->R; r++) {
			printf("X[%d] = %.15e\n", r, mpf_get_d(X[r]));
		}
		
		printf("\nQ (mean queue lengths):\n");
		for (int i = 0; i < qnm->M; i++) {
			printf("Q[%d] =", i+1);
			for (int r = 0; r < qnm->R; r++) {
				printf("\t%.15e", mpf_get_d(Q[i][r]));
			}
			printf("\t(total: %.15e)\n", 
				mpf_get_d(Q[i][0]) + mpf_get_d(Q[i][1]));
		}
		
		printf("=========================================\n");
		printf("Elapsed time (CA): %g s\n", t1-t0);
	}
	mpf_clear(fval);
	
	for (int i=0;i<=qnm->M;i++) mpq_clear(G[i]);
	for (int i=0;i<qnm->R;i++) mpf_clear(X[i]);
	for (int i=0;i<qnm->M;i++) {
		for (int r=0;r<qnm->R;r++) mpf_clear(Q[i][r]);
		free(Q[i]);
	}
	free(G);
	free(X);
	free(Q);
	
	return 0;
}

