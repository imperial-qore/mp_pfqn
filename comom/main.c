#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"

// Define global variables declared in headers
int AORSCTR;
int MULCTR;
int DIVCTR;
double AORSTIME;
double MULTIME;
double DIVTIME;
double t0,t1;
struct rusage ruse;
qnmodel* qnm;

void printcompact(int*n,int R, double elapsed_time)
{
	int s;
	static int first_call = 1;
	if (!first_call) {
		fprintf(stdout,"\r\033[K");
	} else {
		first_call = 0;
	}
	fprintf(stdout,"n=(");
	for (s=0;s<R-1;s++)
	    fprintf(stdout,"%d,",n[s]);
	fprintf(stdout,"%d) - %.6f s",n[R-1], elapsed_time);
	fflush(stdout);
}

int main(int argc, char**argv)
{
	t0=CPUTIME;
	int i,r,nr,s,d;
	int *n;
	combsrep* Dn=NULL;
	combsrep* Dn_old=NULL;
	double step_start, step_elapsed;
	double class_start, class_elapsed;

	/* parse command line arguments */
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	bool debug_output = false;
	int perturbation_digit = 0;
	int perturbation_seed = 23000;
	char* model_file = NULL;
	bool auto_perturbation = false;
	
	if(argc < 2)
	{
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug] [-h|--help] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --debug      : Print performance metrics as exact rational numbers\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		return -1;
	}
	
	for(int i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
			// verbose option not used in COMOM
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
		} else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--debug") == 0) {
			debug_output = true;
		} else if(strcmp(argv[i], "-p") == 0) {
			if(i + 1 < argc) {
				perturbation_digit = atoi(argv[++i]);
				if(perturbation_digit < 1) {
					printf("Error: Perturbation digit must be at least 1\n");
					return -1;
				}
			} else {
				printf("Error: -p option requires a digit argument\n");
				return -1;
			}
		} else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--seed") == 0) {
			if(i + 1 < argc) {
				perturbation_seed = atoi(argv[++i]);
			} else {
				printf("Error: -s option requires a seed argument\n");
				return -1;
			}
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-p digit] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --debug      : Print performance metrics as exact rational numbers\n");
			printf("  -h, --help       : Print this help message\n");
			printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
			return 0;
		} else if(argv[i][0] == '-') {
			// Skip unknown options silently
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug] [-h|--help] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --debug      : Print performance metrics as exact rational numbers\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		return -1;
	}
	
	qnm=(qnmodel*)readmodel(model_file);
	
// solve_attempt: // No longer needed since automatic perturbation is disabled
	// Apply perturbation if requested
	mpz_t scale_factor;
	mpz_init(scale_factor);
	mpz_set_ui(scale_factor, 1);
	
	// Store original values before perturbation for display
	mpz_t** original_L = NULL;
	mpz_t* original_Z = NULL;
	
	if(perturbation_digit > 0) {
		// Store original values before perturbation
		original_L = (mpz_t**)calloc(qnm->M, sizeof(mpz_t*));
		original_Z = (mpz_t*)calloc(qnm->R, sizeof(mpz_t));
		for(int i = 0; i < qnm->M; i++) {
			original_L[i] = (mpz_t*)calloc(qnm->R, sizeof(mpz_t));
			for(int j = 0; j < qnm->R; j++) {
				mpz_init(original_L[i][j]);
				mpz_set(original_L[i][j], qnm->L[i][j]);
			}
		}
		for(int j = 0; j < qnm->R; j++) {
			mpz_init(original_Z[j]);
			mpz_set(original_Z[j], qnm->Z[j]);
		}
		
		// Calculate scale factor: 10^d where d is the perturbation digit
		mpz_set_ui(scale_factor, 10);
		mpz_pow_ui(scale_factor, scale_factor, perturbation_digit);
		
		// For each class, generate a random permutation of 1 to M+1
		int* perm = (int*)calloc(qnm->M + 1, sizeof(int));
		
		for(int j = 0; j < qnm->R; j++) {
			// Generate permutation for class j using seed based on class index
			// Initialize with 1 to M+1
			for(int i = 0; i < qnm->M + 1; i++) {
				perm[i] = i + 1;
			}
			
			// Use a simple linear congruential generator for reproducible randomness
			unsigned int rand_state = perturbation_seed + j * 1000;
			
			// Fisher-Yates shuffle
			for(int i = qnm->M; i > 0; i--) {
				// Generate random number using LCG
				rand_state = rand_state * 1103515245 + 12345;
				int k = (rand_state / 65536) % (i + 1);
				
				// Swap elements i and k
				int temp = perm[i];
				perm[i] = perm[k];
				perm[k] = temp;
			}
			
			// Apply perturbation to L[i][j] for all stations i
			for(int i = 0; i < qnm->M; i++) {
				mpz_mul(qnm->L[i][j], qnm->L[i][j], scale_factor);
				mpz_add_ui(qnm->L[i][j], qnm->L[i][j], perm[i]);
			}
			
			// Apply perturbation to Z[j] using the last element of permutation
			mpz_mul(qnm->Z[j], qnm->Z[j], scale_factor);
			mpz_add_ui(qnm->Z[j], qnm->Z[j], perm[qnm->M]);
		}
		
		free(perm);
		
		if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
			if (auto_perturbation) {
				printf("\nNote: Automatic perturbation applied at digit %d due to singular matrix.\n", perturbation_digit);
			} else {
				printf("\nWarning: Perturbation applied at digit %d. An approximate solution will be computed.\n\n", perturbation_digit);
			}
		}
	}
	
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		if(perturbation_digit > 0) {
			printmodel_with_perturbation(qnm, perturbation_digit, scale_factor, perturbation_seed, original_L, original_Z);
		} else {
			printmodel(qnm);
		}
	}

	/** initializations **/
	nckinit(qnm->M+qnm->R+1,qnm->M); /* declare maximum nck term */

	mpq_vec_t g=NULL;
	mpq_vec_t* g_m=calloc(sizeof(mpq_vec_t),qnm->M+1);
	
	// Storage for marginal normalizing constants G(N-e_r) needed for throughput computation
	mpq_t* marginal_G = (mpq_t*)calloc(qnm->R, sizeof(mpq_t));
	for (int i = 0; i < qnm->R; i++) {
		mpq_init(marginal_G[i]);
	}
	
	// Storage for marginal G^k values needed for queue length computation  
	// Gk_r(N-e_r) for each class r (like in MOM)
	mpq_t* marginal_Gk = (mpq_t*)calloc(qnm->R, sizeof(mpq_t));
	for (int r = 0; r < qnm->R; r++) {
		mpq_init(marginal_Gk[r]);
	}

	/* solve the linear system for all classes */
	for(r=1;r<=qnm->R;r++)
	{
	class_start = CPUTIME;

	/* compute Dn */
	Dn=calloc(sizeof(combsrep),1);
	Dn->n=r;
	Dn->k=qnm->M;
	Dn->card=nck(qnm->M+r-1,qnm->M);
	Dn->combs=(int**)multichoose(Dn->n,Dn->k);
	for (i=1; i<=Dn->card; i++)
		Dn->combs[i-1][r-1]=0;
	Dn->combs=(int**)sortbynnzpos(Dn->combs,Dn->card,Dn->n);

	n=(int*)int_vec(r,0); /* current population vector */
	for (s=1;s<=r-1;s++) 
		n[s-1]=qnm->N[s-1];

	long int cardG = nck(qnm->M+r-1,qnm->M);
	long int cardGk = cardG*qnm->M; 


	g = (mpq_vec_t) mpq_vec(cardG+cardGk,0,1);
	
	int* lu_indices=NULL;
	mpq_vec_t b1 = (mpq_vec_t) mpq_vec(cardGk,0,1);
	mpq_vec_t b1b = (mpq_vec_t) mpq_vec(cardGk,0,1);
	mpq_vec_t G = (mpq_vec_t) mpq_vec(cardG,0,1);

	/* initialize g */
	if (r==1)
	{
		int* comb=(int*)int_vec(r,0); 
		for (i=0;i<=qnm->M;i++)
			mpq_set_si(g[hash(Dn,comb,i+1)-1],1,1);
	}
	else
	{
		for (d=1;d<=Dn->card;d++)
			for (i=0;i<=qnm->M;i++)
			{
				int col=hash(Dn,Dn->combs[d-1],i+1)-1;
				s=Dn->combs[d-1][r-2];
				Dn->combs[d-1][r-2]=0;
				int col_old=hash(Dn_old,Dn->combs[d-1],i+1)-1;
				Dn->combs[d-1][r-2]=s;
				mpq_set(g[col],g_m[s][col_old]);
			}
	}
	LS* linsys=NULL;
		for (nr=1;nr<=qnm->N[r-1];nr++)
		{
			n[r-1] = nr;	
			step_start = CPUTIME;
			step_elapsed = step_start - t0;
			if (!log_output && !normconst_output && !throughput_output && !queue_output) {
				printcompact(n,qnm->R,step_elapsed);
			}
	
			if (n[r-1]==1)
			{
				linsys=setupls(Dn, qnm, n, r);
				lu_indices = mpq_ludcmp(linsys->A11,cardGk);
			}
			/* compute G=B2*g/n */
			mpq_mspvecmul(G, linsys->B2, g);
			mpq_t ninv; mpq_init(ninv); mpq_set_si(ninv,1,n[r-1]); 
			for(i=0;i<cardG;i++) mpq_mul(G[i],ninv,G[i]);

			/* compute b1=B1*g - A12*G */
			mpq_mspvecmul(b1,linsys->B1,g);
			mpq_mspvecmul(b1b,linsys->A12,G);
			for(i=1;i<=cardGk;i++) mpq_sub(b1[i-1],b1[i-1],b1b[i-1]);

			/* compute Gk */
			if (lu_indices == NULL) {
				// LU decomposition indices not computed - matrix is singular
				goto singular_matrix_detected;
			}
			if (mpq_lubksb(linsys->A11, b1, cardGk, lu_indices) < 0) {
				goto singular_matrix_detected;
			}
			copy_Gk_in_g(b1,g);
			copy_G_in_g(G,g);
		
			/* save last M g vectors for initialization of next class */	
			if(r<qnm->R && qnm->N[r-1]-nr<=qnm->M)
			{
				g_m[qnm->N[r-1]-nr]=mpq_vec(cardG+cardGk,0,1);
				mpq_vecdup(g_m[qnm->N[r-1]-nr],g,cardG+cardGk);
				Dn_old = Dn;
			}
			
			/* Save marginal normalizing constant G(N-e_r) for throughput computation */
			/* This is done when nr = N[r-1]-1, i.e., one less job in class r */
			if (nr == qnm->N[r-1] - 1) {
				// G(N-e_r) is stored at position cardGk of the current g vector
				mpq_set(marginal_G[r-1], g[cardGk]);
				
				// Also save G^k value for queue length computation
				// This is equivalent to Gk[r] in MOM - one value per class
				// We need to find the appropriate G^k value from the current state
				// For single class, we can use the first G^k value (queue 1)
				int* comb_reduced = (int*)calloc(r, sizeof(int));
				for (int j = 0; j < r-1; j++) {
					comb_reduced[j] = n[j];
				}
				// Apply the same logic as MOM: for the last class, use the G from before processing it
				if (r == qnm->R) {
					// For the last class, use the previous iteration (nr+1 instead of nr)
					comb_reduced[r-1] = nr + 1; // Use N[r-1] instead of N[r-1]-1
				} else {
					comb_reduced[r-1] = nr; // This is N[r-1]-1
				}
				
				// Get G^k value for the first queue (k=1) which represents Gk[r]
				// In hash function: i=1 means G value, i=2 means G^1 (first queue)
				int gk_idx = hash(Dn, comb_reduced, 2) - 1;
				if (gk_idx >= 0 && gk_idx < cardGk) {
					mpq_set(marginal_Gk[r-1], g[gk_idx]);
				}
				free(comb_reduced);
			}
		}
		class_elapsed = CPUTIME - class_start;
		if (!log_output && !normconst_output && !throughput_output && !queue_output) {
			printf(" (Class %d: %.6f s)\n", r, class_elapsed);
		}	
	}
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\n");
	}
	// Print performance indices
	mpf_t fval;
	mpf_init(fval);
	
	// G(N) is stored at position cardGk in the g vector
	// The first cardGk elements are G^k values, followed by cardG elements which are G values
	long int finalCardGk = nck(qnm->M+qnm->R-1,qnm->M)*qnm->M;
	mpq_t G_total;
	mpq_init(G_total);
	mpq_set(G_total, g[finalCardGk]);
	
	mpf_set_q(fval, G_total);
	double logG = log(mpf_get_d(fval));
	
	// Apply scaling correction to G if perturbation was used
	mpq_t G_scaled;
	mpq_init(G_scaled);
	mpq_set(G_scaled, G_total);
	
	if (mpz_cmp_ui(scale_factor, 1) > 0) {
		// G needs to be divided by scale_factor^Ntot
		int Ntot = 0;
		for (int j = 0; j < qnm->R; j++) {
			Ntot += qnm->N[j];
		}
		mpz_t scale_power;
		mpz_init(scale_power);
		mpz_pow_ui(scale_power, scale_factor, Ntot);
		mpq_t divisor;
		mpq_init(divisor);
		mpq_set_z(divisor, scale_power);
		mpq_div(G_scaled, G_scaled, divisor);
		mpq_clear(divisor);
		mpz_clear(scale_power);
	}
	
	mpf_set_q(fval, G_scaled);
	logG = log(mpf_get_d(fval));
	
	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		// Print exact numerator and denominator
		mpz_t num, den;
		mpz_init(num);
		mpz_init(den);
		mpq_get_num(num, G_scaled);
		mpq_get_den(den, G_scaled);
		gmp_printf("%Zd\n", num);
		gmp_printf("%Zd\n", den);
		mpz_clear(num);
		mpz_clear(den);
	} else if (normconst_g_output) {
		// Print normalizing constant as double
		printf("%.15e\n", mpf_get_d(fval));
	} else if (throughput_output) {
		// Print only throughputs, one per row
		mpq_t X_r;
		mpq_init(X_r);
		for (int r = 1; r <= qnm->R; r++) {
			mpq_div(X_r, marginal_G[r-1], G_total);
			mpq_t X_scaled;
			mpq_init(X_scaled);
			mpq_set(X_scaled, X_r);
			if (mpz_cmp_ui(scale_factor, 1) > 0) {
				mpq_t multiplier;
				mpq_init(multiplier);
				mpq_set_z(multiplier, scale_factor);
				mpq_mul(X_scaled, X_scaled, multiplier);
				mpq_clear(multiplier);
			}
			mpf_set_q(fval, X_scaled);
			printf("%.15e\n", mpf_get_d(fval));
			mpq_clear(X_scaled);
		}
		mpq_clear(X_r);
	} else if (queue_output) {
		// Print queue lengths, all classes for same queue on same row
		mpq_t Q_kr, tmp, tmp2;
		mpq_init(Q_kr);
		mpq_init(tmp);
		mpq_init(tmp2);
		for (int k = 1; k <= qnm->M; k++) {
			for (int r = 1; r <= qnm->R; r++) {
				mpq_set_z(tmp, qnm->L[k-1][r-1]);
				mpq_mul(tmp2, marginal_Gk[r-1], tmp);
				mpq_div(Q_kr, tmp2, G_total);
				mpf_set_q(fval, Q_kr);
				printf("%.15e", mpf_get_d(fval));
				if (r < qnm->R) printf(" ");
			}
			printf("\n");
		}
		mpq_clear(Q_kr);
		mpq_clear(tmp);
		mpq_clear(tmp2);
	} else {
		printf("\n========== Performance Metrics ==========\n");
		
		if (debug_output) {
			// Print exact rational G, omit log(G) in debug mode
			gmp_printf("G = %Qd\n", G_scaled);
		} else {
			// Print double G and log(G)
			printf("G = %.15e\n", mpf_get_d(fval));
			printf("log(G) = %.15e\n", logG);
		}
		
		// Compute exact throughputs using stored marginal normalizing constants
		// X[r] = G(N-e_r) / G(N)
		printf("\nX (throughputs):\n");
		mpq_t X_r;
		mpq_init(X_r);
		
		for (r = 1; r <= qnm->R; r++) {
			// Compute exact throughput: X[r] = G(N-e_r) / G(N)
			mpq_div(X_r, marginal_G[r-1], G_total);
			
			// Apply scaling correction: X[r] needs to be multiplied by scale_factor
			mpq_t X_scaled;
			mpq_init(X_scaled);
			mpq_set(X_scaled, X_r);
			
			if (mpz_cmp_ui(scale_factor, 1) > 0) {
				mpq_t multiplier;
				mpq_init(multiplier);
				mpq_set_z(multiplier, scale_factor);
				mpq_mul(X_scaled, X_scaled, multiplier);
				mpq_clear(multiplier);
			}
			
			if (debug_output) {
				// Print exact rational throughputs
				gmp_printf("X[%d] = %Qd\n", r, X_scaled);
			} else {
				// Print double throughputs
				mpf_set_q(fval, X_scaled);
				printf("X[%d] = %.15e\n", r, mpf_get_d(fval));
			}
			mpq_clear(X_scaled);
		}
		
		
		// Compute queue lengths Q[k,r] = L[k,r] * G^k_r(N-e_r) / G(N)
		printf("\nQ (mean queue lengths):\n");
		
		// Debug: print marginal_Gk values for single class case
		if (qnm->R == 1) {
			printf("Debug: marginal_Gk values:\n");
			mpf_t debug_f; mpf_init(debug_f);
			for (int i = 0; i < qnm->R; i++) {
				mpf_set_q(debug_f, marginal_Gk[i]);
				printf("  marginal_Gk[%d] = %.15e\n", i, mpf_get_d(debug_f));
			}
			mpf_clear(debug_f);
		}
		mpq_t Q_kr, tmp, tmp2, total_q;
		mpq_init(Q_kr);
		mpq_init(tmp);
		mpq_init(tmp2);
		
		for (int k = 1; k <= qnm->M; k++) {
			printf("Q[%d] =", k);
			mpq_init(total_q);
			mpq_set_ui(total_q, 0, 1);
			
			for (r = 1; r <= qnm->R; r++) {
				// Q[k,r] = L[k,r] * Gk[r] / G(N) (same as MOM)
				mpq_set_z(tmp, qnm->L[k-1][r-1]);
				mpq_mul(tmp2, marginal_Gk[r-1], tmp);
				mpq_div(Q_kr, tmp2, G_total);
				mpq_add(total_q, total_q, Q_kr);
				
				if (debug_output) {
					// Print exact rational queue lengths
					printf("\t");
					gmp_printf("%Qd", Q_kr);
				} else {
					// Print double queue lengths
					mpf_set_q(fval, Q_kr);
					printf("\t%.15e", mpf_get_d(fval));
				}
			}
			
			if (debug_output) {
				printf("\t(total: ");
				gmp_printf("%Qd", total_q);
				printf(")\n");
			} else {
				mpf_set_q(fval, total_q);
				printf("\t(total: %.15e)\n", mpf_get_d(fval));
			}
			mpq_clear(total_q);
		}
		
		mpq_clear(Q_kr);
		mpq_clear(tmp);
		mpq_clear(tmp2);
		
		printf("=========================================\n");
		
		mpq_clear(X_r);
	}
	
	mpq_clear(G_total);
	mpq_clear(G_scaled);
	mpf_clear(fval);
	
	// Cleanup marginal_G array
	for (int i = 0; i < qnm->R; i++) {
		mpq_clear(marginal_G[i]);
	}
	free(marginal_G);
	
	// Cleanup marginal_Gk array
	for (int r = 0; r < qnm->R; r++) {
		mpq_clear(marginal_Gk[r]);
	}
	free(marginal_Gk);
	
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		t1=CPUTIME;
		printf("\nElapsed time (CoMoM): %g s\n",t1-t0);
	}

	mpz_clear(scale_factor);
	return 0;

singular_matrix_detected:
	// Report singularity and exit without automatic perturbation
	if (!auto_perturbation && perturbation_digit == 0) {
		fprintf(stderr, "\nError: Model cannot be solved exactly by the CoMoM solver.\n");
		fprintf(stderr, "The system of equations is singular and perturbation is required.\n");
		fprintf(stderr, "Please run manually with perturbation to obtain an approximate solution.\n");
		fprintf(stderr, "\nExample invocations:\n");
		fprintf(stderr, "  %s %s -p 20     # Apply perturbation at 20th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 25     # Apply perturbation at 25th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 15 -s 12345  # Different perturbation at 15th digit with custom seed\n\n", argv[0], model_file);
		mpz_clear(scale_factor);
		return 1;
	} else {
		fprintf(stderr, "\nError: Model is still singular even with perturbation.\n");
		fprintf(stderr, "Please try running manually with a different perturbation digit.\n");
		fprintf(stderr, "\nExample invocations:\n");
		fprintf(stderr, "  %s %s -p 25     # Apply perturbation at 25th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 30     # Apply perturbation at 30th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 15 -s 12345  # Different perturbation at 15th digit with custom seed\n\n", argv[0], model_file);
		mpz_clear(scale_factor);
		return 1;
	}
}

