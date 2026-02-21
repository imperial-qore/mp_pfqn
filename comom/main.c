#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
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
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--debug] [-p digit] [-s seed] model.qn\n", argv[0]);
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
		} else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--debug") == 0) {
			debug_output = true; // -d now prints basis normalizing constants
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
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--debug] [-p digit] [-s seed] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --debug      : Print performance metrics as exact rational numbers\n");
			printf("  -h, --help       : Print this help message\n");
			printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
			printf("  -s seed          : Set perturbation seed (default: 23000)\n");
			return 0;
		} else if(argv[i][0] == '-') {
			// Skip unknown options silently
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--debug] [-p digit] [-s seed] model.qn\n", argv[0]);
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
	
solve_attempt:
	// Apply perturbation if requested
	mpz_t scale_factor;
	mpz_init(scale_factor);
	mpz_set_ui(scale_factor, 1);
	
	// Store original values before perturbation for display
	mpz_t** original_L = NULL;
	mpz_t* original_Z = NULL;
	
	if(perturbation_digit > 0) {
		// Store original values before perturbation
		// Check for potential overflow in allocation
		if (qnm->M < 0 || qnm->R < 0 || 
		    qnm->M > SIZE_MAX / sizeof(mpz_t*) || 
		    qnm->R > SIZE_MAX / sizeof(mpz_t)) {
			fprintf(stderr, "Error: Model dimensions too large for allocation\n");
			exit(EXIT_FAILURE);
		}
		size_t r_size = (size_t)qnm->R;
		size_t m_size = (size_t)qnm->M;
		original_L = (mpz_t**)calloc(m_size, sizeof(mpz_t*));
		original_Z = (mpz_t*)calloc(r_size, sizeof(mpz_t));
		for(int i = 0; i < qnm->M; i++) {
			original_L[i] = (mpz_t*)calloc(r_size, sizeof(mpz_t));
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
	if (g_m == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for g_m array\n");
		return 1;
	}
	
	// Initialize all g_m arrays to prevent accessing uninitialized memory
	// Use the maximum possible size based on the largest class dimensions
	long int max_cardG = nck(qnm->M+qnm->R-1,qnm->M);
	long int max_cardGk = max_cardG * qnm->M;
	long int max_size = max_cardG + max_cardGk;
	
	for (int i = 0; i <= qnm->M; i++) {
		g_m[i] = mpq_vec(max_size, 0, 1);
		if (g_m[i] == NULL) {
			fprintf(stderr, "Error: Memory allocation failed for g_m[%d]\n", i);
			// Clean up previously allocated
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < max_size; k++) {
					mpq_clear(g_m[j][k]);
				}
				free(g_m[j]);
			}
			free(g_m);
			return 1;
		}
	}
	
	// Storage for marginal normalizing constants G(N-e_r) needed for throughput computation
	mpq_t* marginal_G = (mpq_t*)calloc(qnm->R, sizeof(mpq_t));
	if (marginal_G == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for marginal_G array\n");
		free(g_m);
		return 1;
	}
	for (int i = 0; i < qnm->R; i++) {
		mpq_init(marginal_G[i]);
	}
	
	// Storage for marginal G^k values needed for queue length computation  
	// G^k(N-e_r) for each station k and class r - array [M][R]
	mpq_t** marginal_Gk = (mpq_t**)calloc(qnm->M, sizeof(mpq_t*));
	if (marginal_Gk == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for marginal_Gk array\n");
		free(g_m);
		for (int i = 0; i < qnm->R; i++) {
			mpq_clear(marginal_G[i]);
		}
		free(marginal_G);
		return 1;
	}
	for (int k = 0; k < qnm->M; k++) {
		marginal_Gk[k] = (mpq_t*)calloc(qnm->R, sizeof(mpq_t));
		if (marginal_Gk[k] == NULL) {
			fprintf(stderr, "Error: Memory allocation failed for marginal_Gk[%d]\n", k);
			// Clean up
			for (int j = 0; j < k; j++) {
				for (int r = 0; r < qnm->R; r++) {
					mpq_clear(marginal_Gk[j][r]);
				}
				free(marginal_Gk[j]);
			}
			free(marginal_Gk);
			for (int i = 0; i < qnm->R; i++) {
				mpq_clear(marginal_G[i]);
			}
			free(marginal_G);
			free(g_m);
			return 1;
		}
		for (int r = 0; r < qnm->R; r++) {
			mpq_init(marginal_Gk[k][r]);
		}
	}
	
	// Storage for G_1 (previous iteration's G vector) needed for final class performance measures
	mpq_vec_t g_prev = NULL;

	/* solve the linear system for all classes */
	for(r=1;r<=qnm->R;r++)
	{
	class_start = CPUTIME;

	/* compute Dn */
	Dn=calloc(sizeof(combsrep),1);
	if (Dn == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for Dn at class %d\n", r);
		return 1;
	}
	Dn->n=r;
	Dn->k=qnm->M;
	Dn->card=nck(qnm->M+r-1,qnm->M);
	Dn->combs=(int**)multichoose(Dn->n,Dn->k);
	if (Dn->combs == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for Dn->combs at class %d\n", r);
		return 1;
	}
	for (i=1; i<=Dn->card; i++)
		Dn->combs[i-1][r-1]=0;
	Dn->combs=(int**)sortbynnzpos(Dn->combs,Dn->card,Dn->n);

	n=(int*)int_vec(qnm->R,0); /* current population vector */
	for (s=1;s<=r-1;s++) 
		n[s-1]=qnm->N[s-1];

	long int cardG = nck(qnm->M+r-1,qnm->M);
	if (cardG <= 0) {
		fprintf(stderr, "Error: Invalid cardG value %ld at class %d (overflow or calculation error)\n", cardG, r);
		return 1;
	}
	long int cardGk = cardG*qnm->M;
	if (cardGk <= 0 || cardGk < cardG) {
		fprintf(stderr, "Error: Invalid cardGk value %ld at class %d (overflow or calculation error)\n", cardGk, r);
		return 1;
	} 


	g = (mpq_vec_t) mpq_vec(cardG+cardGk,0,1);
	if (g == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for g vector (size %ld) at class %d\n", cardG+cardGk, r);
		return 1;
	}
	
	mpq_vec_t b1 = (mpq_vec_t) mpq_vec(cardGk,0,1);
	if (b1 == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for b1 vector (size %ld) at class %d\n", cardGk, r);
		return 1;
	}
	mpq_vec_t b1b = (mpq_vec_t) mpq_vec(cardGk,0,1);
	if (b1b == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for b1b vector (size %ld) at class %d\n", cardGk, r);
		return 1;
	}
	mpq_vec_t G = (mpq_vec_t) mpq_vec(cardG,0,1);
	if (G == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for G vector (size %ld) at class %d\n", cardG, r);
		return 1;
	}

	/* initialize g */
	if (r==1)
	{
		int* comb=(int*)int_vec(r,0);  // Zero population vector
		// Initialize all elements to 0 first
		for (int idx = 0; idx < cardG + cardGk; idx++) {
			mpq_set_ui(g[idx], 0, 1);
		}
		// Set initial conditions: G(0,m) = 1 for all m
		// m=0 means no queue added (G value)
		// m=1..M means queue m added (G^k values)
		for (i=0;i<=qnm->M;i++) {
			int hash_idx = hash(Dn,comb,i);
			if (hash_idx < 0 || hash_idx > cardG + cardGk) {
				fprintf(stderr, "Error: Invalid hash index %d for i=%d in basis initialization\n", hash_idx, i);
				free(comb);
				return 1;
			}
			mpq_set_si(g[hash_idx-1],1,1);
		}
		free(comb);
	}
	else
	{
		// Create temporary array for hash function with correct size
		int* temp_comb = (int*)calloc(r-1, sizeof(int));
		if (temp_comb == NULL) {
			fprintf(stderr, "Error: Memory allocation failed for temp_comb at class %d\n", r);
			return 1;
		}
		
		for (d=1;d<=Dn->card;d++)
			for (i=0;i<=qnm->M;i++)
			{
				int col=hash(Dn,Dn->combs[d-1],i)-1;
				s=Dn->combs[d-1][r-2];
				
				if (s < 0 || s > qnm->M) {
					fprintf(stderr, "\nError at class %d, d=%d, i=%d:\n", r, d, i);
					fprintf(stderr, "  s = Dn->combs[%d][%d] = %d is out of bounds (should be 0-%d)\n", d-1, r-2, s, qnm->M);
					fprintf(stderr, "  Dn->combs[%d] = [", d-1);
					for (int j = 0; j < r; j++) {
						fprintf(stderr, "%d%s", Dn->combs[d-1][j], j < r-1 ? ", " : "");
					}
					fprintf(stderr, "]\n");
					free(temp_comb);
					return 1;
				}
				
				if (g_m[s] == NULL) {
					fprintf(stderr, "\nError at class %d, d=%d, i=%d:\n", r, d, i);
					fprintf(stderr, "  g_m[%d] is NULL (was never initialized)\n", s);
					fprintf(stderr, "  g_m array status:\n");
					for (int j = 0; j <= qnm->M; j++) {
						fprintf(stderr, "    g_m[%d] = %s\n", j, g_m[j] ? "initialized" : "NULL");
					}
					free(temp_comb);
					return 1;
				}
				
				// Copy first r-1 elements to temp array for hash with Dn_old
				for (int j = 0; j < r-1; j++) {
					temp_comb[j] = Dn->combs[d-1][j];
				}
				temp_comb[r-2] = 0;  // Set the element we're varying to 0
				
				int col_old=hash(Dn_old,temp_comb,i)-1;
				
				// Calculate the old cardG and cardGk for bounds checking
				long int old_cardG = nck(qnm->M+(r-1)-1,qnm->M);
				long int old_cardGk = old_cardG*qnm->M;
				
				if (col_old < 0 || col_old >= old_cardG + old_cardGk) {
					fprintf(stderr, "\nError at class %d, d=%d, i=%d:\n", r, d, i);
					fprintf(stderr, "  col_old = %d is out of bounds (should be 0-%ld)\n", col_old, old_cardG + old_cardGk - 1);
					fprintf(stderr, "  temp_comb = [");
					for (int j = 0; j < r-1; j++) {
						fprintf(stderr, "%d%s", temp_comb[j], j < r-2 ? ", " : "");
					}
					fprintf(stderr, "]\n");
					free(temp_comb);
					return 1;
				}
				
				mpq_set(g[col],g_m[s][col_old]);
			}
		
		free(temp_comb);
	}

	LS* linsys=NULL;
	double setup_start, setup_elapsed;

		/* Save initial g (before any iteration) into g_m if needed.
		 * When N[r] <= M, the inner loop doesn't produce enough saved vectors.
		 * g_m[N[r]] should hold the initial g for class r, corresponding to
		 * population n_r=0.  This matches the MATLAB code's Branch A propagation
		 * equations that carry the initial basis through the full matrix solve. */
		if (r < qnm->R && qnm->N[r-1] <= qnm->M) {
			int idx = qnm->N[r-1];
			mpq_vecdup(g_m[idx], g, cardG+cardGk);
			Dn_old = Dn;
		}

		for (nr=1;nr<=qnm->N[r-1];nr++)
		{
			n[r-1] = nr;	
			step_start = CPUTIME;
			step_elapsed = step_start - t0;
			if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
				printcompact(n,qnm->R,step_elapsed);
			}
	
			if (n[r-1]==1)
			{
				setup_start = CPUTIME;
				int show_progress = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output;
				linsys=setupls(Dn, qnm, n, r, setup_start, show_progress);
				
				// BTF block decomposition and factorization
				linsys->btf = btf_decompose(linsys->A11, Dn, r, qnm->M, cardGk);
				if (!linsys->btf || btf_factorize(linsys->btf) < 0)
					goto singular_matrix_detected;
				setup_elapsed = CPUTIME - setup_start;
				
				if (show_progress) {
					printf("\r\033[K");
					printcompact(n,qnm->R,step_elapsed);
					printf(" [Basis setup: %.6f s]", setup_elapsed);
					fflush(stdout);
				}
			}
			/* compute G=B2*g/n */
			mpq_mspvecmul(G, linsys->B2, g);
			mpq_t ninv; mpq_init(ninv); mpq_set_si(ninv,1,n[r-1]); 
			for(i=0;i<cardG;i++) mpq_mul(G[i],ninv,G[i]);

			/* compute b1=B1*g - A12*G */
			mpq_mspvecmul(b1,linsys->B1,g);
			mpq_mspvecmul(b1b,linsys->A12,G);
			for(i=1;i<=cardGk;i++) mpq_sub(b1[i-1],b1[i-1],b1b[i-1]);

			/* compute Gk via BTF block substitution */
			if (btf_solve(linsys->btf, b1, cardGk) < 0)
				goto singular_matrix_detected;

			/* Fix impossible-combo contamination.
			 *
			 * The A11 matrix couples possible and impossible combo unknowns
			 * through CE shift terms (comb+e_s may exceed N[s]).  The BTF
			 * solve produces a solution where impossible-position values
			 * are non-zero, contaminating possible-position values.
			 *
			 * Fix: build a reduced linear system containing only possible-
			 * combo equations and unknowns, solve it, and replace the
			 * contaminated BTF result. */
			if (r >= 2) {
				/* Identify possible / impossible combos */
				int has_impossible = 0;
				int* combo_possible = (int*)calloc(Dn->card, sizeof(int));
				int n_possible = 0;
				for (d = 1; d <= Dn->card; d++) {
					int is_imp = 0;
					for (s = 0; s <= r-2; s++) {
						if (Dn->combs[d-1][s] > qnm->N[s]) {
							is_imp = 1; break;
						}
					}
					combo_possible[d-1] = !is_imp;
					if (!is_imp) n_possible++;
					if (is_imp) has_impossible = 1;
				}

				if (has_impossible) {
					int M = qnm->M;
					int n_poss_Gk = n_possible * M;

					/* Column mapping: sub-system col → original A11 col */
					int* col_map = (int*)calloc(n_poss_Gk, sizeof(int));
					int ci = 0;
					for (d = 0; d < Dn->card; d++) {
						if (!combo_possible[d]) continue;
						for (int k = 1; k <= M; k++)
							col_map[ci++] = hash(Dn, Dn->combs[d], k) - 1;
					}

					/* Zero G at impossible positions */
					for (d = 0; d < Dn->card; d++)
						if (!combo_possible[d])
							mpq_set_ui(G[d], 0, 1);

					/* Recompute RHS = B1*g - A12*G  (G zeroed at impossible) */
					mpq_vec_t b1_rhs = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
					mpq_vec_t a12g   = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
					mpq_mspvecmul(b1_rhs, linsys->B1, g);
					mpq_mspvecmul(a12g,   linsys->A12, G);
					for (i = 0; i < cardGk; i++)
						mpq_sub(b1_rhs[i], b1_rhs[i], a12g[i]);

					/* Count included rows.  Replicate setupls iteration
					 * order to identify which original A11 row corresponds
					 * to each CE/PC equation.
					 *   Include CE  rows for POSSIBLE combos.
					 *   Include PC  rows for POSSIBLE combos where the
					 *               shifted combo d+e_s is also POSSIBLE.  */
					int* row_included = (int*)calloc(cardGk, sizeof(int));
					int n_sub_rows = 0;
					{
						int row = 0;
						for (d = 1; d <= Dn->card; d++) {
							if (int_vecsubsum(Dn->combs[d-1], 0, r-1) >= M)
								continue;
							int dp = combo_possible[d-1];
							/* M CE rows */
							for (int k = 0; k < M; k++) {
								if (dp) { row_included[row] = 1; n_sub_rows++; }
								row++;
							}
							/* r-1 PC rows */
							for (s = 0; s < r-1; s++) {
								if (dp && Dn->combs[d-1][s] + 1 <= qnm->N[s]) {
									row_included[row] = 1;
									n_sub_rows++;
								}
								row++;
							}
						}
					}

					/* Build sub-matrix  A_sub (n_sub_rows × n_poss_Gk)
					 * and sub-RHS       b_sub (n_sub_rows)              */
					mpq_mat_t A_sub = (mpq_mat_t)mpq_matzeros(n_sub_rows, n_poss_Gk);
					mpq_vec_t b_sub = (mpq_vec_t) mpq_vec(n_sub_rows, 0, 1);
					{
						int sr = 0;
						for (int orig_row = 0; orig_row < cardGk; orig_row++) {
							if (!row_included[orig_row]) continue;
							for (int j = 0; j < n_poss_Gk; j++)
								mpq_set(A_sub[sr][j],
								        linsys->A11[orig_row][col_map[j]]);
							mpq_set(b_sub[sr], b1_rhs[orig_row]);
							sr++;
						}
					}

					/* Gaussian elimination with partial pivoting
					 * on the  m × n  system  (m ≥ n).                */
					for (int col = 0; col < n_poss_Gk; col++) {
						int piv = -1;
						for (int row = col; row < n_sub_rows; row++) {
							if (mpq_sgn(A_sub[row][col]) != 0)
								{ piv = row; break; }
						}
						if (piv < 0) {
							fprintf(stderr,
							  "Warning: singular possible sub-system at col %d\n", col);
							break;
						}
						if (piv != col) {
							for (int j = 0; j < n_poss_Gk; j++)
								mpq_swap(A_sub[col][j], A_sub[piv][j]);
							mpq_swap(b_sub[col], b_sub[piv]);
						}
						for (int row = col + 1; row < n_sub_rows; row++) {
							if (mpq_sgn(A_sub[row][col]) == 0) continue;
							mpq_t fac; mpq_init(fac);
							mpq_div(fac, A_sub[row][col], A_sub[col][col]);
							for (int j = col; j < n_poss_Gk; j++) {
								mpq_t tmp; mpq_init(tmp);
								mpq_mul(tmp, fac, A_sub[col][j]);
								mpq_sub(A_sub[row][j], A_sub[row][j], tmp);
								mpq_clear(tmp);
							}
							mpq_t tmp; mpq_init(tmp);
							mpq_mul(tmp, fac, b_sub[col]);
							mpq_sub(b_sub[row], b_sub[row], tmp);
							mpq_clear(tmp);
							mpq_clear(fac);
						}
					}

					/* Back-substitution */
					mpq_vec_t x_sub = (mpq_vec_t) mpq_vec(n_poss_Gk, 0, 1);
					for (int col = n_poss_Gk - 1; col >= 0; col--) {
						mpq_set(x_sub[col], b_sub[col]);
						for (int j = col + 1; j < n_poss_Gk; j++) {
							mpq_t tmp; mpq_init(tmp);
							mpq_mul(tmp, A_sub[col][j], x_sub[j]);
							mpq_sub(x_sub[col], x_sub[col], tmp);
							mpq_clear(tmp);
						}
						mpq_div(x_sub[col], x_sub[col], A_sub[col][col]);
					}

					/* Copy solution into b1 for possible combos */
					for (int j = 0; j < n_poss_Gk; j++)
						mpq_set(b1[col_map[j]], x_sub[j]);

					/* Zero b1 at impossible combo positions */
					for (d = 0; d < Dn->card; d++)
						if (!combo_possible[d])
							for (int k = 1; k <= M; k++)
								mpq_set_ui(b1[hash(Dn, Dn->combs[d], k) - 1], 0, 1);

					/* Cleanup */
					for (i = 0; i < n_sub_rows; i++) {
						for (int j = 0; j < n_poss_Gk; j++)
							mpq_clear(A_sub[i][j]);
						free(A_sub[i]);
					}
					free(A_sub);
					for (i = 0; i < n_sub_rows; i++) mpq_clear(b_sub[i]);
					free(b_sub);
					for (i = 0; i < n_poss_Gk; i++) mpq_clear(x_sub[i]);
					free(x_sub);
					for (i = 0; i < cardGk; i++) { mpq_clear(b1_rhs[i]); mpq_clear(a12g[i]); }
					free(b1_rhs); free(a12g);
					free(col_map);
					free(row_included);
				}
				free(combo_possible);
			}

			/* Save G_1 for final class - MATLAB: G_1=G before G=Gn */
			/* This happens at the LAST iteration of the LAST class */
			if (r == qnm->R && nr == qnm->N[r-1]) {
				// Save current g before overwriting with new values
				if (g_prev == NULL) {
					g_prev = (mpq_vec_t) mpq_vec(cardG+cardGk,0,1);
				}
				if (g_prev != NULL) {
					mpq_vecdup(g_prev, g, cardG+cardGk);
				}
			}

			copy_Gk_in_g(b1,g);
			copy_G_in_g(G,g);

				/* save last M g vectors for initialization of next class */
			if(r<qnm->R && qnm->N[r-1]-nr<=qnm->M)
			{
				int idx = qnm->N[r-1]-nr;
				if (idx < 0 || idx > qnm->M) {
					fprintf(stderr, "Error: Invalid g_m index %d at class %d, nr=%d\n", idx, r, nr);
					return 1;
				}
				// Copy g to g_m[idx] - no need to allocate since already initialized
				mpq_vecdup(g_m[idx],g,cardG+cardGk);
				Dn_old = Dn;
			}
			
			/* Save marginal normalizing constants for performance measures */
			if (r == qnm->R && nr == qnm->N[r-1] - 1) {
				mpq_set(marginal_G[qnm->R-1], g[cardGk]); // G(N-e_R) for class R
				// Save G^k(N-e_R) for each station k
				for (int k = 0; k < qnm->M; k++) {
					mpq_set(marginal_Gk[k][qnm->R-1], g[k]);
				}
			} else if (r == qnm->R && nr == qnm->N[r-1]) {
				int d = 1;
				int * comb = calloc(sizeof(int),r);

				for (int s = 1; s <= r; s++)
					comb[s-1]=Dn->combs[d-1][s-1];

				for (int s = 1; s <= qnm->R-1; s++) {
					comb[s-1]++;

					// Calculate index using hash function - 0 for G values
					int index = hash(Dn, comb, 0) - 1;

					mpq_set(marginal_G[s-1], g[index]);

					// Save G^k values for this class
					for (int k = 1; k <= qnm->M; k++) {
						index = hash(Dn, comb, k) - 1;
						// Print all contents of marginal_G and marginal_Gk
						mpq_set(marginal_Gk[k-1][s-1], g[index]);
					}
					comb[s-1]--;
				}
				free(comb);

				// For final class R, save penultimate basis for debug output
				// Free existing allocation
				//for (int t=0; t<cardG+cardGk; t++) mpq_clear(g_m[1][t]);
				//free(g_m[1]);
				//g_m[1] = mpq_vec(cardG+cardGk,0,1);
				//mpq_vecdup(g_m[1], g, cardG+cardGk);
			}
		}
		class_elapsed = CPUTIME - class_start;

		if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
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
		mpq_canonicalize(G_scaled);
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
			// Use the same logic as in normal output mode
			if (r < qnm->R) {
				// For classes s=1 to R-1: use marginal_G[r-1]
				mpq_div(X_r, marginal_G[r-1], G_total);
			} else {
				// For final class R: use g_prev if available
				if (g_prev != NULL) {
					long int finalCardGk = nck(qnm->M + qnm->R - 1, qnm->M) * qnm->M;
					mpq_div(X_r, g_prev[finalCardGk], G_total);
				} else {
					// Fallback if g_prev not available
					mpq_div(X_r, marginal_G[r-1], G_total);
				}
			}
			
			// Apply scaling correction if perturbation was used
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
				// Check if the original demand was 0.0
				mpz_t original_L_value;
				mpz_init(original_L_value);
				mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
				
				if (mpz_cmp_ui(original_L_value, 0) == 0) {
					// Original demand was 0, so Q should be 0
					mpq_set_ui(Q_kr, 0, 1);
				} else {
					// Q[k,r] = L[k,r] * G^k(N-e_r) / G(N)
					mpq_set_z(tmp, qnm->L[k-1][r-1]);
					
					if (r < qnm->R) {
						// For classes s=1 to R-1: use the saved marginal_Gk values
						mpq_mul(tmp2, marginal_Gk[k-1][r-1], tmp);
						mpq_div(Q_kr, tmp2, G_total);
					} else {
						// For final class R: use g_prev if available
						if (g_prev != NULL) {
							// G^k(N) for station k is at position (k-1) in g_prev
							int gk_position = (k-1);
							mpq_mul(tmp2, g_prev[gk_position], tmp);
							mpq_div(Q_kr, tmp2, G_total);
						} else {
							// Fallback if g_prev not available
							mpq_mul(tmp2, marginal_Gk[k-1][r-1], tmp);
							mpq_div(Q_kr, tmp2, G_total);
						}
					}
				}
				
				mpz_clear(original_L_value);
				
				// Multiply by the multiplicity mi[k-1] for station k
				if (qnm->mi[k-1] > 1) {
					mpq_t mi_q;
					mpq_init(mi_q);
					mpq_set_ui(mi_q, qnm->mi[k-1], 1);
					mpq_mul(Q_kr, Q_kr, mi_q);
					mpq_clear(mi_q);
				}
				
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
		printf("========== Performance Metrics ==========\n");
		
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
		// For classes 1 to R-1, compute G(N-e_r) from the final basis
		// We need to search the final g vector for the right normalizing constants
		// The final basis has all combinations with total population N
		
		// For classes 1 to R-1, we need G(N-e_r) from the final basis
		// The challenge is that in multi-class models, the g vector contains
		// normalizing constants for different queue population distributions,
		// not directly for class populations.
		// 
		// The combinations are sorted by number of zeros (more zeros first),
		// then by leftmost zero position. We need to find the right combination
		// that represents N-e_r for each class.
		
		printf("\nX (throughputs):\n");
		mpq_t X_r;
		mpq_init(X_r);
		
		for (r = 1; r <= qnm->R; r++) {
			// Compute exact throughput using MATLAB approach
			if (r < qnm->R) {
				// For classes s=1 to R-1: X(s) = G(hash(N,oner(N,s),0+1))/ G(hash(N,N,0+1))
				// 
				// TODO: For classes 1 to R-1, need to extract G(N-e_r) from final basis
				// Currently using marginal_G which was saved at wrong population for multi-class
				mpq_div(X_r, marginal_G[r-1], G_total);
			} else {
				// For final class R: X(R) = G_1(hash(N,N,0+1))/ G(hash(N,N,0+1))
				// G_1 is the g vector from iteration N[R]-1, use its G(N) value
				if (g_prev != NULL) {
					long int finalCardGk = nck(qnm->M + qnm->R - 1, qnm->M) * qnm->M;
					// G(N) is at position finalCardGk in g_prev
					
					
					mpq_div(X_r, g_prev[finalCardGk], G_total);
				} else {
					// Fallback if g_prev not available
					mpq_div(X_r, marginal_G[r-1], G_total);
				}
			}
			
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
		
		mpq_t Q_kr, tmp, tmp2, total_q;
		mpq_init(Q_kr);
		mpq_init(tmp);
		mpq_init(tmp2);
		
		for (int k = 1; k <= qnm->M; k++) {
			printf("Q[%d] =", k);
			mpq_init(total_q);
			mpq_set_ui(total_q, 0, 1);
			
			for (r = 1; r <= qnm->R; r++) {
				// Check if the original demand was 0.0
				mpz_t original_L_value;
				mpz_init(original_L_value);
				mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
				
				if (mpz_cmp_ui(original_L_value, 0) == 0) {
					// Original demand was 0, so Q should be 0
					mpq_set_ui(Q_kr, 0, 1);
				} else {
					// Compute queue lengths using MATLAB approach
					mpq_set_z(tmp, qnm->L[k-1][r-1]);
					
					if (r < qnm->R) {
						// For classes s=1 to R-1: Q(k,s) = L(k,s)*G(hash(N,oner(N,s),k+1))/ G(hash(N,N,0+1))
						// Use the saved marginal_Gk values for station k and class r
						mpq_mul(tmp2, marginal_Gk[k-1][r-1], tmp);
						mpq_div(Q_kr, tmp2, G_total);
					} else {
						// For final class R: Q(k,R) = L(k,R)*G_1(hash(N,N,k+1))/ G(hash(N,N,0+1))
						// Use g_prev (G_1) at position corresponding to G^k(N)
						if (g_prev != NULL) {
							// G^k(N) for station k is at position (k-1) in g_prev
							// The first M elements of g_prev are G^1(N), G^2(N), ..., G^M(N)
							int gk_position = (k-1);
							mpq_mul(tmp2, g_prev[gk_position], tmp);
							mpq_div(Q_kr, tmp2, G_total);
						} else {
							// Fallback if g_prev not available
							mpq_mul(tmp2, marginal_Gk[k-1][r-1], tmp);
							mpq_div(Q_kr, tmp2, G_total);
						}
					}
				}
				
				mpz_clear(original_L_value);
				
				// Multiply by the multiplicity mi[k-1] for station k
				if (qnm->mi[k-1] > 1) {
					mpq_t mi_q;
					mpq_init(mi_q);
					mpq_set_ui(mi_q, qnm->mi[k-1], 1);
					mpq_mul(Q_kr, Q_kr, mi_q);
					mpq_clear(mi_q);
				}
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
		t1=CPUTIME;
		printf("Elapsed time (CoMoM): %g s\n",t1-t0);
		
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
	
	// Cleanup marginal_Gk 2D array [M][R]
	for (int k = 0; k < qnm->M; k++) {
		for (int r = 0; r < qnm->R; r++) {
			mpq_clear(marginal_Gk[k][r]);
		}
		free(marginal_Gk[k]);
	}
	free(marginal_Gk);
	
	// Cleanup g_m arrays
	for (int i = 0; i <= qnm->M; i++) {
		if (g_m[i]) {
			for (int k = 0; k < max_size; k++) {
				mpq_clear(g_m[i][k]);
			}
			free(g_m[i]);
		}
	}
	free(g_m);
	
	mpz_clear(scale_factor);
	return 0;

singular_matrix_detected:
	// If singular and no perturbation yet, apply automatic perturbation at digit 20
	if (!auto_perturbation && perturbation_digit == 0) {
		// If exact normalizing constant was requested, don't auto-perturb
		if (normconst_output) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by the CoMoM solver.\n");
			fprintf(stderr, "The system of equations is singular. No exact normalizing constant available.\n");
			mpz_clear(scale_factor);
			return 1;
		}
		fprintf(stderr, "\nWarning: Model cannot be solved exactly by the CoMoM solver.\n");
		fprintf(stderr, "The system of equations is singular. Automatically applying perturbation at digit 20.\n");
		perturbation_digit = 20;
		auto_perturbation = true;
		
		// Store original parameters if not already stored
		if (original_L == NULL) {
			// Check for potential overflow in allocation
			if (qnm->M < 0 || qnm->R < 0 || 
			    qnm->M > SIZE_MAX / sizeof(mpz_t*) || 
			    qnm->R > SIZE_MAX / sizeof(mpz_t)) {
				fprintf(stderr, "Error: Model dimensions too large for allocation\n");
				exit(EXIT_FAILURE);
			}
			size_t r_size = (size_t)qnm->R;
			size_t m_size = (size_t)qnm->M;
			original_L = (mpz_t**)calloc(m_size, sizeof(mpz_t*));
			original_Z = (mpz_t*)calloc(r_size, sizeof(mpz_t));
			for(int i = 0; i < qnm->M; i++) {
				original_L[i] = (mpz_t*)calloc(r_size, sizeof(mpz_t));
				for(int j = 0; j < qnm->R; j++) {
					mpz_init(original_L[i][j]);
					mpz_set(original_L[i][j], qnm->L[i][j]);
				}
			}
			for(int j = 0; j < qnm->R; j++) {
				mpz_init(original_Z[j]);
				mpz_set(original_Z[j], qnm->Z[j]);
			}
		}
		
		// Re-read the model to reset to original values
		qnm=(qnmodel*)readmodel(model_file);
		
		// Clear and reset variables for retry
		mpz_clear(scale_factor);
		// Free allocated memory for g_m and marginal arrays before retry
		if (g_m) {
			for(int i = 0; i <= qnm->M; i++) {
				if(g_m[i]) {
					// Clear all mpq_t values in g_m[i] before freeing
					for (int k = 0; k < max_size; k++) {
						mpq_clear(g_m[i][k]);
					}
					free(g_m[i]);
				}
			}
			free(g_m);
		}
		if (marginal_G) {
			for (int i = 0; i < qnm->R; i++) {
				mpq_clear(marginal_G[i]);
			}
			free(marginal_G);
		}
		if (marginal_Gk) {
			for (int k = 0; k < qnm->M; k++) {
				if (marginal_Gk[k]) {
					for (int r = 0; r < qnm->R; r++) {
						mpq_clear(marginal_Gk[k][r]);
					}
					free(marginal_Gk[k]);
				}
			}
			free(marginal_Gk);
		}
		
		// Jump back to retry with perturbation
		goto solve_attempt;
		
	} else if (auto_perturbation) {
		// Automatic perturbation already tried and still singular
		fprintf(stderr, "\nError: Model is still singular even with automatic perturbation at digit 20.\n");
		fprintf(stderr, "Please try running manually with a different perturbation digit.\n");
		fprintf(stderr, "\nExample invocations:\n");
		fprintf(stderr, "  %s %s -p 25     # Apply perturbation at 25th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 30     # Apply perturbation at 30th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 15 -s 12345  # Different perturbation at 15th digit with custom seed\n\n", argv[0], model_file);
		mpz_clear(scale_factor);
		return 1;
	} else {
		// Manual perturbation was used but still singular
		fprintf(stderr, "\nError: Model is still singular even with perturbation at digit %d.\n", perturbation_digit);
		fprintf(stderr, "Please try running manually with a different perturbation digit.\n");
		fprintf(stderr, "\nExample invocations:\n");
		fprintf(stderr, "  %s %s -p 25     # Apply perturbation at 25th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 30     # Apply perturbation at 30th digit\n", argv[0], model_file);
		fprintf(stderr, "  %s %s -p 15 -s 12345  # Different perturbation at 15th digit with custom seed\n\n", argv[0], model_file);
		mpz_clear(scale_factor);
		return 1;
	}
}

