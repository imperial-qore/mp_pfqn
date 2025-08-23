#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"

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
	DEBUG=false;
	VERBOSE=true;
	int t,i,r,s,Ntot,h,k;
	int *n;
	int cardG, cardGk;
	mpq_t tmp; 
	mpq_t tmp2; 
	mpq_vec_t g=NULL,Gk=NULL,G=NULL,b=NULL,b2=NULL,gr=NULL;
	LS* A=NULL;
	double step_start, step_elapsed;
	double class_start, class_elapsed;

	/* parse command line arguments */
	bool verbose_output = false;
	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool throughput_output = false;
	bool queue_output = false;
	bool debug_output = false;
	bool bounds_output = false;
	int perturbation_digit = 0;
	int perturbation_seed = 23000;
	char* model_file = NULL;
	
	if(argc < 2)
	{
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug] [-b|--bounds] [-h|--help] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --debug      : Enable debug output (progress messages, timing, etc.)\n");
		printf("  -b, --bounds     : Compute performance bounds (requires -p, shows ranges)\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		return -1;
	}
	
	for(int i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
			verbose_output = true;
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
		} else if(strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bounds") == 0) {
			bounds_output = true;
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug] [-b|--bounds] [-h|--help] [-p digit] [-s seed] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --debug      : Enable debug output (progress messages, timing, etc.)\n");
			printf("  -h, --help       : Print this help message\n");
			printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
			printf("  -s seed          : Set perturbation seed (default: 23000)\n");
			return 0;
		} else if(strcmp(argv[i], "-p") == 0) {
			if(i + 1 < argc) {
				perturbation_digit = atoi(argv[++i]);
				if(perturbation_digit < 1 || perturbation_digit > 15) {
					printf("Error: Perturbation digit must be between 1 and 15\n");
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
		} else if(argv[i][0] == '-') {
			// Skip unknown options silently - no additional logic needed
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug] [-b|--bounds] [-h|--help] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --debug      : Enable debug output (progress messages, timing, etc.)\n");
		printf("  -b, --bounds     : Compute performance bounds (requires -p, shows ranges)\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		return -1;
	}
	
	// Validate bounds option
	if(bounds_output && perturbation_digit == 0) {
		printf("Error: -b/--bounds option requires -p option to be specified\n");
		return -1;
	}
	
	qnm=(qnmodel*)readmodel(model_file);
	
	// Store original parameters for printing
	int** original_L = NULL;
	int* original_Z = NULL;
	if(perturbation_digit > 0) {
		// Allocate memory for original parameters
		original_L = (int**)malloc(qnm->M * sizeof(int*));
		for(int i = 0; i < qnm->M; i++) {
			original_L[i] = (int*)malloc(qnm->R * sizeof(int));
			for(int j = 0; j < qnm->R; j++) {
				original_L[i][j] = qnm->L[i][j];
			}
		}
		original_Z = (int*)malloc(qnm->R * sizeof(int));
		for(int j = 0; j < qnm->R; j++) {
			original_Z[j] = qnm->Z[j];
		}
	}
	
	// Apply perturbation if requested
	long scale_factor = 1;
	if(perturbation_digit > 0) {
		// Calculate scale factor: 10^d where d is the perturbation digit
		for(int i = 0; i < perturbation_digit; i++) {
			scale_factor *= 10;
		}
		
		// Scale all parameters
		for(int i = 0; i < qnm->M; i++) {
			for(int j = 0; j < qnm->R; j++) {
				// Add small perturbation using improved randomization
				long hash_val = perturbation_seed * 1103515245 + i * 12345 + j * 67891;
				long perturb = (abs(hash_val) % 9) + 1; // 1-9 range
				qnm->L[i][j] = qnm->L[i][j] * scale_factor + perturb;
			}
		}
		
		for(int j = 0; j < qnm->R; j++) {
			qnm->Z[j] = qnm->Z[j] * scale_factor;
		}
		
		if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
			if (bounds_output) {
				printf("\nWarning: Perturbation applied at digit %d. Bounds will be computed using ±perturbation.\n\n", perturbation_digit);
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
	/* initializations */
	nckinit(qnm->M+qnm->R-1,qnm->R); /* declare maximum nck term */
	n=(int*)int_vec(qnm->R,0); /* current population vector */
	mpq_init(tmp);
	mpq_init(tmp2);
	int m=qnm->M;
	Ntot=0;
	/* solve the linear system for all classes */
	for(r=1;r<=qnm->R;r++)
	{
	class_start = CPUTIME;
	Ntot += qnm->N[r-1];
	cardGk=nck(m+r-1,r)*r;
	cardG=nck(m+r-2,r-1)*r;
	
	if(r==1) /* if this is the first processed population */
	{
		/* initialize linear system */	
		A =(LS*) setupls(qnm->L,n,qnm->Z,qnm->mi,m,r);
		if (A == NULL) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by the MoM solver.\n");
			fprintf(stderr, "The system of equations is singular. Consider introducing small perturbations\n");
			fprintf(stderr, "to the service demands to allow approximate solution.\n");
			fprintf(stderr, "\nExample invocations:\n");
			fprintf(stderr, "  %s %s -p 6     # Apply perturbation at 6th digit\n", argv[0], model_file);
			fprintf(stderr, "  %s %s -p 12    # Apply perturbation at 12th digit\n\n", argv[0], model_file);
			return 1;
		}
		b=(mpq_vec_t)mpq_vec(m,0,1);
		b2=(mpq_vec_t)mpq_vec(m,0,1);
		/* the norm.consts for n=(0,0,...,0) are equal to 1 */
		g=(mpq_vec_t)mpq_vec(m+1,1,1); 
		gr=(mpq_vec_t)mpq_vec(m+1,1,1);
		Gk=(mpq_vec_t)mpq_vec(m,1,1);
		G=(mpq_vec_t)mpq_vec(1,1,1);
	}
	else
	{
		t=0;
		for (h=1;h<=A->H;h++)
		for (k=0;k<nck(m,h);k++)
		{
			free(A->C[t].diag);
			free(A->C[t++].nondiag->coeff);
		}
		free(A->A12); free(A->B1r); free(A->B2r); free(A);
		A =(LS*) setupls(qnm->L,n,qnm->Z,qnm->mi,m,r);
		if (A == NULL) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by the MoM solver.\n");
			fprintf(stderr, "The system of equations is singular. Consider introducing small perturbations\n");
			fprintf(stderr, "to the service demands to allow approximate solution.\n");
			fprintf(stderr, "\nExample invocations:\n");
			fprintf(stderr, "  %s %s -p 6     # Apply perturbation at 6th digit\n", argv[0], model_file);
			fprintf(stderr, "  %s %s -p 12    # Apply perturbation at 12th digit\n\n", argv[0], model_file);
			return 1;
		}
		for (t=0;t<nck(m+r-2,r-1)*(r-1);t++)
		{
			mpq_clear(b[t]);
			mpq_clear(b2[t]);
		}
		for (t=0;t<nck(m+r-3,r-2)*(r-1);t++)
			mpq_clear(G[t]);
		b=(mpq_vec_t)mpq_vec(cardGk,0,1);
		b2=(mpq_vec_t)mpq_vec(cardGk,0,1);
		G=(mpq_vec_t)mpq_vec(cardG,0,1);
		Gk=(mpq_vec_t)mpq_vec(cardGk,0,1);
		for (i=0;i<nck(m+r-2,r-1);i++)
		{
			for (s=0;s<=r-2;s++)
			{
				mpq_set(G[i*r+s],g[i*(r-1)+s]);
			}
			mpq_set(G[i*r+r-1],gr[i*(r-1)]);
		}
		mpq_mspvecmul(b,A->A12,G);
		for(t=0;t<cardGk;t++)
			mpq_neg(b[t],b[t]);	
		Gk=(mpq_vec_t)blocksolve(A,b);
		if (Gk == NULL) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by the MoM solver.\n");
			fprintf(stderr, "The system of equations is singular. Consider introducing small perturbations\n");
			fprintf(stderr, "to the service demands to allow approximate solution.\n");
			fprintf(stderr, "\nExample invocations:\n");
			fprintf(stderr, "  %s %s -p 6     # Apply perturbation at 6th digit\n", argv[0], model_file);
			fprintf(stderr, "  %s %s -p 12    # Apply perturbation at 12th digit\n\n", argv[0], model_file);
			return 1;
		}
		for (t=0;t<nck(m+r-2,r-1)*(r-1);t++)
		{
			mpq_clear(g[t]);
			mpq_clear(gr[t]);
		}
		g=(mpq_vec_t)mpq_vec(cardGk+cardG,0,1);
		gr=(mpq_vec_t)mpq_vec(cardGk+cardG,0,1);
		copy_Gk_in_g(Gk,g);
		copy_G_in_g(G,g);
		/* annihilate unused normalizing constants */
		for(s=0;s<r;s++) if(qnm->Z[s]==0) for(t=0;t<nck(m+r-2,r-1);t++) mpq_set_si(g[cardGk+t*r+1+s],0,1);
	}
	/* start solving the sequence of linear systems */
	for (n[r-1]=1;n[r-1]<=qnm->N[r-1];n[r-1]++) /* for all population of class r */ 
	{
		step_start = CPUTIME;
		step_elapsed = step_start - t0;
		if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
			printcompact(n,qnm->R,step_elapsed);
		}
		/** at this point g is equal to gr **/
		/* if this is the last iteration, save gr */
		if(n[r-1]==qnm->N[r-1]) mpq_vecdup(gr,g,cardGk+cardG);

		/* compute G=B2r*gr/nr */

		mpq_mspvecmul(G,A->B2r,g);
		mpq_t ninv; mpq_init(ninv); mpq_set_si(ninv,1,n[r-1]); 
		for(t=0;t<cardG;t++) mpq_mul(G[t],ninv,G[t]);

		/* compute b=B1r*gr - A12*G */
		mpq_mspvecmul(b,A->B1r,g);
		mpq_mspvecmul(b2,A->A12,G);
		for(t=0;t<cardGk;t++) mpq_sub(b[t],b[t],b2[t]);

		/* compute Gk */
		Gk=(mpq_vec_t)blocksolve(A,b);
		if (Gk == NULL) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by then MoM solver.\n");
			fprintf(stderr, "The system of equations is singular. Consider introducing small perturbations\n");
			fprintf(stderr, "to the service demands to allow approximate solution.\n");
			fprintf(stderr, "\nExample invocations:\n");
			fprintf(stderr, "  %s %s -p 6     # Apply perturbation at 6th digit\n", argv[0], model_file);
			fprintf(stderr, "  %s %s -p 12    # Apply perturbation at 12th digit\n\n", argv[0], model_file);
			return 1;
		}
		
		/** after this point g is no longer equal to gr **/
		copy_Gk_in_g(Gk,g);
		copy_G_in_g(G,g);
		free_Gk();
		/* annihilate unused normalizing constants */
		for(s=0;s<r;s++) if(qnm->Z[s]==0) for(t=0;t<nck(m+r-2,r-1);t++) mpq_set_si(g[cardGk+t*r+1+s],0,1);
		/* G,g and gr are not freed because we need them at the next cycle */	
	}
	n[r-1]=qnm->N[r-1];
	step_elapsed = CPUTIME - t0;
	if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printcompact(n,qnm->R,step_elapsed);
		class_elapsed = CPUTIME - class_start;
		printf(" (Class %d: %.6f s)\n", r, class_elapsed);
	}
	}
	if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\n");
	}
	mdecrease(qnm,G,Gk,g,gr,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
	mpq_clear(tmp);
	mpq_clear(tmp2);
	t1=CPUTIME;
	if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\nElapsed time (MOM): %.6f s\n",t1-t0);
	}

	// Clean up original parameter storage
	if(perturbation_digit > 0) {
		for(int i = 0; i < qnm->M; i++) {
			free(original_L[i]);
		}
		free(original_L);
		free(original_Z);
	}

	return 0;
}

