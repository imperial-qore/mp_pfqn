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

void recal_multi_exact(qnmodel* qn, mpq_t G);
void recal_marginal_exact(qnmodel* qn, int class_to_reduce, mpq_t G_marginal);
void recal_queue_marginal_exact(qnmodel* qn, int class_to_reduce, int station_k, mpq_t G_k_marginal);

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
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
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
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
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
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
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
		// Compute and print throughputs using marginal normalizing constants
		mpq_t G_marginal, X_r;
		mpf_t X_float;
		mpq_init(G_marginal);
		mpq_init(X_r);
		mpf_init(X_float);
		
		for (int r = 0; r < qn->R; r++) {
			// Compute G(N-e_r) for class r
			recal_marginal_exact(qn, r, G_marginal);
			
			// Compute X[r] = G(N-e_r) / G(N)
			mpq_div(X_r, G_marginal, G_exact);
			
			// Convert to double and print
			mpf_set_q(X_float, X_r);
			printf("%.15e\n", mpf_get_d(X_float));
		}
		
		mpq_clear(G_marginal);
		mpq_clear(X_r);
		mpf_clear(X_float);
	} else if (queue_output) {
		// Compute and print queue lengths using marginal normalizing constants
		// Q[k,r] = m[k] * L[k,r] * G^k(N-e_r) / G(N)
		mpq_t G_k_marginal, Q_kr, L_kr, mi_k;
		mpf_t Q_float;
		mpq_init(G_k_marginal);
		mpq_init(Q_kr);
		mpq_init(L_kr);
		mpq_init(mi_k);
		mpf_init(Q_float);
		
		for (int k = 0; k < qn->M; k++) {
			for (int r = 0; r < qn->R; r++) {
				// Compute G^k(N-e_r) for station k and class r
				recal_queue_marginal_exact(qn, r, k, G_k_marginal);
				
				// Q[k,r] = m[k] * L[k,r] * G^k(N-e_r) / G(N)
				mpq_set_ui(mi_k, qn->mi[k], 1);
				mpq_set_z(L_kr, qn->L[k][r]);
				mpq_mul(Q_kr, mi_k, L_kr);
				mpq_mul(Q_kr, Q_kr, G_k_marginal);
				mpq_div(Q_kr, Q_kr, G_exact);
				
				// Convert to double and print
				mpf_set_q(Q_float, Q_kr);
				printf("%.15e", mpf_get_d(Q_float));
				if (r < qn->R - 1) printf(" ");
			}
			printf("\n");
		}
		
		mpq_clear(G_k_marginal);
		mpq_clear(Q_kr);
		mpq_clear(L_kr);
		mpq_clear(mi_k);
		mpf_clear(Q_float);
	} else {
		printf("\n========== Performance Metrics ==========\n");
		printf("Normalizing constant:\n");
		printf("G = %.15e\n", G);
		printf("log(G) = %.15e\n", logG);
		
		// Compute and print throughputs using marginal normalizing constants
		printf("\nX (throughputs):\n");
		mpq_t G_marginal, X_r;
		mpf_t X_float;
		mpq_init(G_marginal);
		mpq_init(X_r);
		mpf_init(X_float);
		
		for (int r = 0; r < qn->R; r++) {
			// Compute G(N-e_r) for class r
			recal_marginal_exact(qn, r, G_marginal);
			
			// Compute X[r] = G(N-e_r) / G(N)
			mpq_div(X_r, G_marginal, G_exact);
			
			// Convert to double and print
			mpf_set_q(X_float, X_r);
			printf("X[%d] = %.15e\n", r+1, mpf_get_d(X_float));
		}
		
		mpq_clear(G_marginal);
		mpq_clear(X_r);
		mpf_clear(X_float);
		
		// Compute and print queue lengths using marginal normalizing constants
		printf("\nQ (mean queue lengths):\n");
		mpq_t G_k_marginal, Q_kr, L_kr, mi_k;
		mpf_t Q_float;
		mpq_init(G_k_marginal);
		mpq_init(Q_kr);
		mpq_init(L_kr);
		mpq_init(mi_k);
		mpf_init(Q_float);
		
		for (int k = 0; k < qn->M; k++) {
			printf("Q[%d] =", k+1);
			mpq_t total_q; 
			mpq_init(total_q); 
			mpq_set_ui(total_q, 0, 1);
			
			for (int r = 0; r < qn->R; r++) {
				// Compute G^k(N-e_r) for station k and class r
				recal_queue_marginal_exact(qn, r, k, G_k_marginal);
				
				// Q[k,r] = m[k] * L[k,r] * G^k(N-e_r) / G(N)
				mpq_set_ui(mi_k, qn->mi[k], 1);
				mpq_set_z(L_kr, qn->L[k][r]);
				mpq_mul(Q_kr, mi_k, L_kr);
				mpq_mul(Q_kr, Q_kr, G_k_marginal);
				mpq_div(Q_kr, Q_kr, G_exact);
				
				// Convert to double and print
				mpf_set_q(Q_float, Q_kr);
				printf("\t%.15e", mpf_get_d(Q_float));
				mpq_add(total_q, total_q, Q_kr);
			}
			mpf_t tval; 
			mpf_init(tval); 
			mpf_set_q(tval, total_q);
			printf("\t(total: %.15e)\n", mpf_get_d(tval));
			mpf_clear(tval);
			mpq_clear(total_q);
		}
		
		mpq_clear(G_k_marginal);
		mpq_clear(Q_kr);
		mpq_clear(L_kr);
		mpq_clear(mi_k);
		mpf_clear(Q_float);
		
		printf("=========================================\n");
	}
	
	mpq_clear(G_exact);
	
	t1 = CPUTIME;
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\nElapsed time (RECAL): %g s\n", t1-t0);
	}
	
	return 0;
}


