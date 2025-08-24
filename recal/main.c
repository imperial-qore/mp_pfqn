#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

void recal_multi(qnmodel* qn, double *X);
void recal_multi_exact(qnmodel* qn, mpq_t G);

int main(int argc,char ** argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;
	
	/* parse command line arguments */
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	char* model_file = NULL;
	
	if(argc < 2)
	{
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log      : Print only log of normalizing constant as double\n");
		printf("  -e, --ex       : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}
	
	for(int i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
			// verbose option not used in RECAL
		} else if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--log") == 0) {
			log_output = true;
		} else if(strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ex") == 0) {
			normconst_output = true;
		} else if(strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--nc") == 0) {
			normconst_g_output = true;
		} else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tput") == 0) {
			throughput_output = true;
		} else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qlen") == 0) {
			queue_output = true;
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log      : Print only log of normalizing constant as double\n");
			printf("  -e, --ex       : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -h, --help       : Print this help message\n");
			return 0;
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log      : Print only log of normalizing constant as double\n");
		printf("  -e, --ex       : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}
	
	qnmodel* qn = (qnmodel*)readmodel(model_file);
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printmodel(qn);
	}

	mpq_t G_exact;
	mpq_init(G_exact);
	recal_multi_exact(qn, G_exact);
	
	// Check if recal failed (G_exact would be zero)
	if (mpq_cmp_ui(G_exact, 0, 1) == 0) {
		mpq_clear(G_exact);
		return -1;
	}
	
	double G = mpq_get_d(G_exact);
	double logG = log(G);
	
	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		// Print exact numerator and denominator
		mpz_t num, den;
		mpz_init(num);
		mpz_init(den);
		mpq_get_num(num, G_exact);
		mpq_get_den(den, G_exact);
		gmp_printf("%Zd\n", num);
		gmp_printf("%Zd\n", den);
		mpz_clear(num);
		mpz_clear(den);
	} else if (normconst_g_output) {
		// Print normalizing constant as double
		printf("%.15e\n", G);
	} else if (throughput_output) {
		// RECAL doesn't compute throughputs, print zeros
		int r;
		for (r=0; r<qn->R; r++) {
			printf("0.000000000000000e+00\n");
		}
	} else if (queue_output) {
		// RECAL doesn't compute queue lengths, print zeros
		int m;
		for (m=0; m<qn->M; m++) {
			int r;
			for (r=0; r<qn->R; r++) {
				printf("0.000000000000000e+00");
				if (r < qn->R - 1) printf(" ");
			}
			printf("\n");
		}
	} else {
		printf("\n========== Performance Metrics ==========\n");
		printf("G = %.15e\n", G);
		printf("log(G) = %.15e\n", logG);
		printf("=========================================\n");
	}
	
	mpq_clear(G_exact);
	
	t1 = CPUTIME;
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\nElapsed time (RECAL): %g s\n", t1-t0);
	}
	
	return 0;
}


