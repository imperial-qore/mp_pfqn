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
	bool exact_output = false;
	char* model_file = NULL;
	
	if(argc < 2)
	{
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
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
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
			printf("  -h, --help       : Print this help message\n");
			return 0;
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--exact] [-h|--help] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -h, --help       : Print this help message\n");
		return -1;
	}
	
	qnmodel* qn = (qnmodel*)readmodel(model_file);
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !exact_output) {
		printmodel(qn);
	}
	
	// Check total population and warn if too large
	int total_population = 0;
	for (int r = 0; r < qn->R; r++) {
		total_population += qn->N[r];
	}
	if (total_population > 30) {
		fprintf(stderr, "\nWarning: Total population is %d (> 30). RECAL is likely to fail due to memory constraints.\n", total_population);
		fprintf(stderr, "Consider using a different solver for models with large populations.\n\n");
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
		mpq_canonicalize(G_exact);
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
	} else if (exact_output) {
		// Print all performance metrics in exact rational form
		printf("========== Performance Metrics (Exact) ==========\n");
		printf("Normalizing constant:\n");
		gmp_printf("G = %Qd\n", G_exact);
		printf("log(G) = %.15e\n", logG);
		
		// Compute and print throughputs using marginal normalizing constants
		printf("\nX (throughputs):\n");
		mpq_t G_marginal, X_r;
		mpq_init(G_marginal);
		mpq_init(X_r);
		
		for (int r = 0; r < qn->R; r++) {
			// Compute G(N-e_r) for class r
			recal_marginal_exact(qn, r, G_marginal);
			
			// Compute X[r] = G(N-e_r) / G(N)
			mpq_div(X_r, G_marginal, G_exact);
			
			// Print exact rational throughput
			gmp_printf("X[%d] = %Qd\n", r+1, X_r);
		}
		
		mpq_clear(G_marginal);
		mpq_clear(X_r);
		
		// Compute and print queue lengths using marginal normalizing constants
		printf("\nQ (mean queue lengths):\n");
		mpq_t G_k_marginal, Q_kr, L_kr, mi_k;
		mpq_init(G_k_marginal);
		mpq_init(Q_kr);
		mpq_init(L_kr);
		mpq_init(mi_k);
		
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
				
				// Print exact rational queue length
				gmp_printf("\t%Qd", Q_kr);
				mpq_add(total_q, total_q, Q_kr);
			}
			gmp_printf("\t(total: %Qd)\n", total_q);
			mpq_clear(total_q);
		}
		
		mpq_clear(G_k_marginal);
		mpq_clear(Q_kr);
		mpq_clear(L_kr);
		mpq_clear(mi_k);
		
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (RECAL): %.6f s\n", t1-t0);
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
		// First compute all metrics
		
		// Allocate arrays to store computed throughputs
		double* X_values = malloc(qn->R * sizeof(double));
		mpq_t G_marginal, X_r;
		mpf_t X_float;
		mpq_init(G_marginal);
		mpq_init(X_r);
		mpf_init(X_float);
		
		// Compute throughputs
		for (int r = 0; r < qn->R; r++) {
			// Compute G(N-e_r) for class r
			recal_marginal_exact(qn, r, G_marginal);
			
			// Compute X[r] = G(N-e_r) / G(N)
			mpq_div(X_r, G_marginal, G_exact);
			
			// Convert to double and store
			mpf_set_q(X_float, X_r);
			X_values[r] = mpf_get_d(X_float);
		}
		
		mpq_clear(G_marginal);
		mpq_clear(X_r);
		mpf_clear(X_float);
		
		// Allocate arrays to store computed queue lengths
		double** Q_values = malloc(qn->M * sizeof(double*));
		double* Q_totals = malloc(qn->M * sizeof(double));
		for (int k = 0; k < qn->M; k++) {
			Q_values[k] = malloc(qn->R * sizeof(double));
		}
		
		// Compute queue lengths
		mpq_t G_k_marginal, Q_kr, L_kr, mi_k;
		mpf_t Q_float;
		mpq_init(G_k_marginal);
		mpq_init(Q_kr);
		mpq_init(L_kr);
		mpq_init(mi_k);
		mpf_init(Q_float);
		
		for (int k = 0; k < qn->M; k++) {
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
				
				// Convert to double and store
				mpf_set_q(Q_float, Q_kr);
				Q_values[k][r] = mpf_get_d(Q_float);
				mpq_add(total_q, total_q, Q_kr);
			}
			mpf_t tval; 
			mpf_init(tval); 
			mpf_set_q(tval, total_q);
			Q_totals[k] = mpf_get_d(tval);
			mpf_clear(tval);
			mpq_clear(total_q);
		}
		
		mpq_clear(G_k_marginal);
		mpq_clear(Q_kr);
		mpq_clear(L_kr);
		mpq_clear(mi_k);
		mpf_clear(Q_float);
		
		// Clear progress line before printing results
		fprintf(stderr, "\r                                                                                  \r");
		// Now print all results
		printf("========== Performance Metrics ==========\n");
		printf("Normalizing constant:\n");
		printf("G = %.15e\n", G);
		printf("log(G) = %.15e\n", logG);
		
		// Print throughputs
		printf("\nX (throughputs):\n");
		for (int r = 0; r < qn->R; r++) {
			printf("X[%d] = %.15e\n", r+1, X_values[r]);
		}
		
		// Print queue lengths
		printf("\nQ (mean queue lengths):\n");
		for (int k = 0; k < qn->M; k++) {
			printf("Q[%d] =", k+1);
			for (int r = 0; r < qn->R; r++) {
				printf("\t%.15e", Q_values[k][r]);
			}
			printf("\t(total: %.15e)\n", Q_totals[k]);
		}
		
		printf("=========================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (RECAL): %g s\n", t1-t0);
		
		// Free allocated memory
		free(X_values);
		for (int k = 0; k < qn->M; k++) {
			free(Q_values[k]);
		}
		free(Q_values);
		free(Q_totals);
	}
	
	mpq_clear(G_exact);
	
	return 0;
}


