#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"
#include "lvtab.h"
#include "singular.h"
#include "../mommod/mmengine.h"

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
bool USE_LINBOX;
mpq_t* MOM_CAPTURE = NULL;


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

void printcompact_with_setup(int*n,int R, double elapsed_time, double setup_time)
{
	int s;
	fprintf(stdout,"\r\033[K");
	fprintf(stdout,"n=(");
	for (s=0;s<R-1;s++)
	    fprintf(stdout,"%d,",n[s]);
	fprintf(stdout,"%d) - %.6f s [Basis setup: %.6f s]",n[R-1], elapsed_time, setup_time);
	fflush(stdout);
}

// Generate a random permutation of integers 1 to n using Fisher-Yates shuffle
void generate_permutation(int* perm, int n, unsigned int seed) {
	// Initialize with 1 to n
	for(int i = 0; i < n; i++) {
		perm[i] = i + 1;
	}
	
	// Use a simple linear congruential generator for reproducible randomness
	unsigned int rand_state = seed;
	
	// Fisher-Yates shuffle
	for(int i = n - 1; i > 0; i--) {
		// Generate random number using LCG
		rand_state = rand_state * 1103515245 + 12345;
		int j = (rand_state / 65536) % (i + 1);
		
		// Swap elements i and j
		int temp = perm[i];
		perm[i] = perm[j];
		perm[j] = temp;
	}
}

void apply_perturbation_to_model(qnmodel* qnm, int perturbation_digit, int perturbation_seed, mpz_t scale_factor) {
	// Calculate scale factor: 10^d where d is the perturbation digit
	mpz_set_ui(scale_factor, 10);
	mpz_pow_ui(scale_factor, scale_factor, perturbation_digit);
	
	// For each class, generate a random permutation of 1 to M+1
	int* perm = (int*)calloc(qnm->M + 1, sizeof(int));
	
	for(int j = 0; j < qnm->R; j++) {
		// Generate permutation for class j using seed based on class index
		generate_permutation(perm, qnm->M + 1, perturbation_seed + j * 1000);
		
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
}

/* Lower replica levels 0..r-2 of the full MoM basis: see util/lvtab.h. */
static void lv_free(mpq_vec_t v, long n)
{
	long t;
	if (!v) return;
	for (t=0;t<=n;t++) mpq_clear(v[t]);
	free(v);
}

/* Advance levels 0..r-2 from population n-1_r to n by the class-r
 * population constraint.  glow points at the old level r-1 (second part
 * of g), which must not yet be overwritten. */
static void lower_step(qnmodel* qnm, lvtab* T, int r, mpq_vec_t* low, mpq_t* glow, long nr)
{
	int l, i, j, s, m = T->m;
	mpq_t acc, c, w, zr;
	mpq_init(acc); mpq_init(c); mpq_init(w); mpq_init(zr);
	mpq_set_z(zr, qnm->Z[r-1]);
	for (l=0;l<=r-2;l++) {
		mpq_t* above = (l+1 == r-1) ? glow : low[l+1];
		for (i=0;i<T->card[l];i++) {
			for (s=0;s<r;s++) {
				mpq_mul(acc, zr, low[l][i*r+s]);
				for (j=0;j<m;j++) {
					if (mpz_sgn(qnm->L[j][r-1]) == 0) continue;
					mpq_set_z(c, qnm->L[j][r-1]);
					mpq_set_si(w, qnm->mi[j] + T->cb[l][i][j], 1);
					mpq_mul(c, c, w);
					mpq_mul(w, c, above[T->up[l][i*m+j]*r+s]);
					mpq_add(acc, acc, w);
				}
				mpq_set_si(w, nr, 1);
				mpq_div(low[l][i*r+s], acc, w);
			}
		}
	}
	mpq_clear(acc); mpq_clear(c); mpq_clear(w); mpq_clear(zr);
}

int solve_model(qnmodel* qnm, int argc, char** argv, bool verbose_output, bool log_output, 
                bool normconst_output, bool normconst_g_output, bool throughput_output, 
                bool queue_output, bool debug_output, bool bounds_output, 
                int perturbation_digit, int perturbation_seed, bool auto_perturbation,
                mpz_t scale_factor, mpz_t** original_L, mpz_t* original_Z) {
	int t,i,r,s,Ntot,h,k;
	int *n;
	int cardG, cardGk;
	mpq_t tmp; 
	mpq_t tmp2; 
	mpq_vec_t g=NULL,Gk=NULL,G=NULL,b=NULL,b2=NULL,gr=NULL;
	LS* A=NULL;
	double step_start, step_elapsed;
	double class_start, class_elapsed;
	double setup_start, setup_elapsed;
	
	mpq_init(tmp);
	mpq_init(tmp2);
	
	/* initializations */
	nckinit(qnm->M+qnm->R-1,qnm->R); /* declare maximum nck term */
	n=(int*)int_vec(qnm->R,0); /* current population vector */
	int m=qnm->M;
	Ntot=0;
	int R=qnm->R, l;
	int full = mom_mixed_think(qnm); /* else the mdecrease descent is exact */
	lvtab* T = full ? lvtab_new(m,R) : NULL;
	mpq_vec_t *low  = (mpq_vec_t*)calloc(R>1?R-1:1,sizeof(mpq_vec_t)); /* levels 0..r-2 */
	mpq_vec_t *lowr = (mpq_vec_t*)calloc(R>1?R-1:1,sizeof(mpq_vec_t)); /* their snapshot at N-1_r */
	
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
			setup_start = CPUTIME;
			int show_progress = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !bounds_output;
			A =(LS*) setupls_progress(qnm->L,n,qnm->Z,qnm->mi,m,r,setup_start,show_progress);
			if (A == NULL) {
				return -1; // Signal singular matrix
			}
			setup_elapsed = CPUTIME - setup_start;
			b=(mpq_vec_t)mpq_vec(m,0,1);
			b2=(mpq_vec_t)mpq_vec(m,0,1);
			/* the norm.consts for n=(0,0,...,0) are equal to 1 */
			g=(mpq_vec_t)mpq_vec(m+1,1,1); 
			gr=(mpq_vec_t)mpq_vec(m+1,0,1); /* G(n-1_1)=0 at n=0: stays 0 if N_1=0 */
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
			setup_start = CPUTIME;
			int show_progress = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !bounds_output;
			A =(LS*) setupls_progress(qnm->L,n,qnm->Z,qnm->mi,m,r,setup_start,show_progress);
			if (A == NULL) {
				return -1; // Signal singular matrix
			}
			setup_elapsed = CPUTIME - setup_start;
			/* seed levels 0..r-2 for class r: entries for classes
			 * 1..r-2 carry over, the class r-1 entry is G(n-1_{r-1})
			 * taken from the snapshot at N-1_{r-1} */
			{
				int cardGk_old = nck(m+r-2,r-1)*(r-1);
				for (l=0;full && l<=r-2;l++) {
					mpq_t* old  = (l == r-2) ? g  + cardGk_old : low[l];
					mpq_t* snap = (l == r-2) ? gr + cardGk_old : lowr[l];
					mpq_vec_t nw = (mpq_vec_t)mpq_vec(T->card[l]*r,0,1);
					for (i=0;i<T->card[l];i++) {
						for (s=0;s<=r-2;s++) mpq_set(nw[i*r+s], old[i*(r-1)+s]);
						mpq_set(nw[i*r+r-1], snap[i*(r-1)]);
					}
					if (l < r-2) { lv_free(low[l],T->card[l]*(r-1)); lv_free(lowr[l],T->card[l]*(r-1)); }
					low[l]  = nw;
					lowr[l] = (mpq_vec_t)mpq_vec(T->card[l]*r,0,1); /* stays 0 if N_r=0 */
				}
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
				return -1; // Signal singular matrix
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
		}
		/* start solving the sequence of linear systems */
		for (n[r-1]=1;n[r-1]<=qnm->N[r-1];n[r-1]++) /* for all population of class r */ 
		{
			step_start = CPUTIME;
			step_elapsed = step_start - t0;
			if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !bounds_output) {
				printcompact_with_setup(n,qnm->R,step_elapsed,setup_elapsed);
			}
			/** at this point g is equal to gr **/
			/* if this is the last iteration, save gr */
			if(n[r-1]==qnm->N[r-1]) {
				mpq_vecdup(gr,g,cardGk+cardG);
				for (l=0;full && l<=r-2;l++) mpq_vecdup(lowr[l],low[l],T->card[l]*r);
			}

			/* compute G=B2r*gr/nr */

			mpq_mspvecmul(G,A->B2r,g);
			mpq_t ninv; mpq_init(ninv); mpq_set_si(ninv,1,n[r-1]); 
			for(t=0;t<cardG;t++) mpq_mul(G[t],ninv,G[t]);
			if (full) lower_step(qnm,T,r,low,g+cardGk,n[r-1]);

			/* compute b=B1r*gr - A12*G */
			mpq_mspvecmul(b,A->B1r,g);
			mpq_mspvecmul(b2,A->A12,G);
			for(t=0;t<cardGk;t++) mpq_sub(b[t],b[t],b2[t]);

			/* compute Gk */
			Gk=(mpq_vec_t)blocksolve(A,b);
			if (Gk == NULL) {
				return -1; // Signal singular matrix
			}
			
			/** after this point g is no longer equal to gr **/
			copy_Gk_in_g(Gk,g);
			copy_G_in_g(G,g);
			free_Gk();
			/* G,g and gr are not freed because we need them at the next cycle */	
		}
		n[r-1]=qnm->N[r-1];
		step_elapsed = CPUTIME - t0;
		if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !bounds_output) {
			printcompact_with_setup(n,qnm->R,step_elapsed,setup_elapsed);
			class_elapsed = CPUTIME - class_start;
			printf(" (Class %d: %.6f s)\n", r, class_elapsed);
		}
	}
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output && !bounds_output) {
		printf("\n");
	}
	if (!full) mdecrease(qnm,G,Gk,g,gr,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
	else
	/* read levels 0 and 1 at N into the layout perfindices expects:
	 * G = [G(N), G(N-1_1..R)], Gk[k] = [G^{+k}(N), G^{+k}(N-1_1..R)] */
	{
		mpq_t *L0, *S0, *L1, *S1;
		int cardGkR = nck(m+R-1,R)*R;
		L0 = (R >= 2) ? low[0]  : g + cardGkR;
		S0 = (R >= 2) ? lowr[0] : gr + cardGkR;
		L1 = (R >= 3) ? low[1]  : (R == 2 ? g + cardGkR  : g);
		S1 = (R >= 3) ? lowr[1] : (R == 2 ? gr + cardGkR : gr);
		mpq_vec_t G0 = (mpq_vec_t)mpq_vec(R+1,0,1);
		mpq_vec_t G1 = (mpq_vec_t)mpq_vec(m*(R+1),0,1);
		for (s=0;s<R;s++) mpq_set(G0[s], L0[s]);
		mpq_set(G0[R], S0[0]);
		for (k=0;k<m;k++) {
			for (s=0;s<R;s++) mpq_set(G1[k*(R+1)+s], L1[k*R+s]);
			mpq_set(G1[k*(R+1)+R], S1[k*R]);
		}
		perfindices(qnm,G0,G1,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
		lv_free(G0,R+1); lv_free(G1,m*(R+1));
	}
	if (full) for (l=0;l<R-1;l++) { lv_free(low[l],T->card[l]*R); lv_free(lowr[l],T->card[l]*R); }
	free(low); free(lowr);
	lvtab_free(T);
	mpq_clear(tmp);
	mpq_clear(tmp2);
	
	return 0; // Success
}

/* Exact solve of a derived model for the tiered fallback (util/singular.h):
 * the constants are captured in out[] instead of printed.  0 on success,
 * 1 if the MoM system of q is singular. */
static int sg_direct(qnmodel* q, mpq_t* out, void* ctx)
{
	mpz_t one; int rc;
	(void)ctx;
	mpz_init_set_ui(one, 1);
	MOM_CAPTURE = out;
	rc = solve_model(q, 0, NULL, false, false, false, false, false, true, false, false,
	                 0, 0, false, one, NULL, NULL);
	MOM_CAPTURE = NULL;
	mpz_clear(one);
	return rc == 0 ? 0 : 1;
}

/* Print the reconstructed constants through perfindices, in the layout it
 * expects: G = [G(N), G(N-1_1..R)], Gk[k] = [G^{+k}(N), G^{+k}(N-1_1..R)]
 * (G^{+k}(N) is not read). */
static void report_exact(qnmodel* q, mpq_t* ex, bool verbose_output, bool log_output,
                         bool normconst_output, bool normconst_g_output, bool throughput_output,
                         bool queue_output, bool debug_output, mpz_t scale_factor)
{
	int R = q->R, M = q->M, i, k;
	mpq_vec_t G0 = (mpq_vec_t)mpq_vec(R+1,0,1);
	mpq_vec_t G1 = (mpq_vec_t)mpq_vec(M*(R+1),0,1);
	for (i=0;i<=R;i++) mpq_set(G0[i], ex[i]);
	for (k=0;k<M;k++) for (i=0;i<R;i++) mpq_set(G1[k*(R+1)+1+i], ex[SG_IDX_GK(R,k,i)]);
	perfindices(q,G0,G1,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,false,scale_factor);
	lv_free(G0,R+1); lv_free(G1,M*(R+1));
}

int main(int argc, char**argv)
{
	t0=CPUTIME;
	DEBUG=false;
	VERBOSE=true;

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
	bool auto_perturbation = false;
	bool use_linbox = false;
	/* -X: never fall back to the automatic perturbation.  On a singular
	 * system the solver reports the singularity and exits non-zero instead
	 * of returning an approximate answer.  This is the fast-fail probe used
	 * by bin/safe_mom to decide when to move to the next tier; an explicit
	 * -p is still honoured. */
	bool no_perturb = false;
	int mom_threads = 1;
	/* --no-reduce: skip the exact reduction (identical stations folded,
	 * proportional classes aggregated) done before every solve */
	bool no_reduce = false;

	if(argc < 2)
	{
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--exact] [-b|--bounds] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -b, --bounds     : Compute performance bounds (requires -p, shows ranges)\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		printf("  --linbox         : Use LinBox Dixon p-adic lifting solver (avoids coefficient blowup)\n");
		printf("  -X, --no-perturb : Fail on a singular system instead of auto-perturbing (exact-only probe)\n");
		printf("  -j n             : Threads for the modular interpolation of singular models (default 1)\n");
		printf("  --no-reduce      : Do not fold identical stations or merge proportional classes first\n");
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
			debug_output = true; // -d enables exact output in perfindices
		} else if(strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bounds") == 0) {
			bounds_output = true;
		} else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-v|--verbose] [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--exact] [-b|--bounds] [-p digit] [-s seed] model.qn\n", argv[0]);
			printf("  -v, --verbose    : Print exact ratios for all performance measures\n");
			printf("  -l, --log        : Print only log of normalizing constant as double\n");
			printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
			printf("  -g, --nc         : Print normalizing constant as double\n");
			printf("  -t, --tput       : Print only throughputs, one per row\n");
			printf("  -q, --qlen       : Print only queue lengths, one per row\n");
			printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
			printf("  -h, --help       : Print this help message\n");
			printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
			printf("  -s seed          : Set perturbation seed (default: 23000)\n");
			printf("  --linbox         : Use LinBox Dixon p-adic lifting solver (avoids coefficient blowup)\n");
			printf("  -X, --no-perturb : Fail on a singular system instead of auto-perturbing (exact-only probe)\n");
			printf("  -j n             : Threads for the modular interpolation of singular models (default 1)\n");
			printf("  --no-reduce      : Do not fold identical stations or merge proportional classes first\n");
			return 0;
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
		} else if(strcmp(argv[i], "-j") == 0) {
			if(i + 1 < argc) mom_threads = atoi(argv[++i]);
		} else if(strcmp(argv[i], "--linbox") == 0) {
			use_linbox = true;
		} else if(strcmp(argv[i], "-X") == 0 || strcmp(argv[i], "--no-perturb") == 0) {
			no_perturb = true;
		} else if(strcmp(argv[i], "--no-reduce") == 0) {
			no_reduce = true;
		} else if(argv[i][0] == '-') {
			// Skip unknown options silently - no additional logic needed
		} else {
			model_file = argv[i];
		}
	}
	
	if(model_file == NULL) {
		printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-h|--help] [-d|--exact] [-b|--bounds] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  -l, --log        : Print only log of normalizing constant as double\n");
		printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
		printf("  -g, --nc         : Print normalizing constant as double\n");
		printf("  -t, --tput       : Print only throughputs, one per row\n");
		printf("  -q, --qlen       : Print only queue lengths, one per row\n");
		printf("  -d, --exact      : Print all performance metrics in full exact precision (integer or rational)\n");
		printf("  -b, --bounds     : Compute performance bounds (requires -p, shows ranges)\n");
		printf("  -h, --help       : Print this help message\n");
		printf("  -p digit         : Apply perturbation at the specified digit (e.g., -p 5)\n");
		printf("  -s seed          : Set perturbation seed (default: 23000)\n");
		printf("  --linbox         : Use LinBox Dixon p-adic lifting solver (avoids coefficient blowup)\n");
		printf("  -X, --no-perturb : Fail on a singular system instead of auto-perturbing (exact-only probe)\n");
		printf("  -j n             : Threads for the modular interpolation of singular models (default 1)\n");
		printf("  --no-reduce      : Do not fold identical stations or merge proportional classes first\n");
		return -1;
	}
	
	// Validate bounds option
	if(bounds_output && perturbation_digit == 0) {
		printf("Error: -b/--bounds option requires -p option to be specified\n");
		return -1;
	}
	
	USE_LINBOX = use_linbox;
	qnm=(qnmodel*)readmodel(model_file);
	qn_require_closed(qnm, "mom", 0);
	
	// Store original parameters for printing
	mpz_t** original_L = NULL;
	mpz_t* original_Z = NULL;
	if(perturbation_digit > 0) {
		// Allocate memory for original parameters
		original_L = (mpz_t**)malloc(qnm->M * sizeof(mpz_t*));
		for(int i = 0; i < qnm->M; i++) {
			original_L[i] = (mpz_t*)malloc(qnm->R * sizeof(mpz_t));
			for(int j = 0; j < qnm->R; j++) {
				mpz_init(original_L[i][j]);
				mpz_set(original_L[i][j], qnm->L[i][j]);
			}
		}
		original_Z = (mpz_t*)malloc(qnm->R * sizeof(mpz_t));
		for(int j = 0; j < qnm->R; j++) {
			mpz_init(original_Z[j]);
			mpz_set(original_Z[j], qnm->Z[j]);
		}
	}
	
	// Apply perturbation if requested
	mpz_t scale_factor;
	mpz_init(scale_factor);
	mpz_set_ui(scale_factor, 1);  // Default to 1 when no scaling
	if(perturbation_digit > 0) {
		apply_perturbation_to_model(qnm, perturbation_digit, perturbation_seed, scale_factor);
		if (normconst_output || debug_output)
			fprintf(stderr, "Warning: -p %d perturbs the model; the rational printed is exact for the perturbed model only and is NOT the exact G of the input model (it changes with -s).\n", perturbation_digit);
		
		if (debug_output && !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
			if (auto_perturbation) {
				printf("\nNote: Automatic perturbation applied at digit %d due to singular matrix.\n", perturbation_digit);
			} else if (bounds_output) {
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
	
	/* Exact reduction first (util/singular.h): stations with identical
	 * demand rows fold into one of summed multiplicity and classes with
	 * proportional demands and think times aggregate, so the MoM basis is
	 * built over fewer stations and classes; the constants of the input
	 * model are recovered exactly.  Not under -p or -b, which perturb. */
	if (!no_reduce && perturbation_digit == 0 && !bounds_output) {
		bool report = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output;
		int nres = SG_NRES(qnm);
		mpq_t* exact = sg_mpq_arr(nres);
		char tier[256];
		if (!no_perturb) mommod_config(mom_threads, 61, 0, 0, 64);
		sg_ops ops = { sg_direct, mommod_interp, NULL, " per prime" };
		int rc = sg_solve_reduced(qnm, exact, tier, sizeof tier, &ops, !no_perturb);
		if (rc != -1) {
			if (rc != 0) {
				if (no_perturb)
					fprintf(stderr, "\nError: Model cannot be solved exactly by the MoM solver.\n"
					                "The system of equations of the reduced model is singular; -X suppresses the singular fallback.\n");
				else
					fprintf(stderr, "\nError: the MoM system of the reduced model is singular, and every exact fallback failed.\n"
					                "bin/ca solves any well-posed model by convolution.\n");
				sg_mpq_arr_free(exact, nres);
				return 1;
			}
			fprintf(stderr, "mom: solved the reduced model (%s)\n", tier);
			if (report) printf("\nSolved the exactly reduced model (%s)\n", tier);
			report_exact(qnm, exact, verbose_output, log_output, normconst_output, normconst_g_output,
			             throughput_output, queue_output, debug_output, scale_factor);
			sg_mpq_arr_free(exact, nres);
			mpz_clear(scale_factor);
			return 0;
		}
		sg_mpq_arr_free(exact, nres);
	}

	// Try to solve the model
	int result = solve_model(qnm, argc, argv, verbose_output, log_output,
	                        normconst_output, normconst_g_output, throughput_output,
	                        queue_output, debug_output, bounds_output,
	                        perturbation_digit, perturbation_seed, auto_perturbation,
	                        scale_factor, original_L, original_Z);

	// If singular and no perturbation yet, apply automatic perturbation
	if (result == -1 && !auto_perturbation && perturbation_digit == 0) {
		// -X: exact-only probe for bin/safe_mom, report the singularity and
		// fail fast instead of running the tiered fallback below
		if (no_perturb) {
			fprintf(stderr, "\nError: Model cannot be solved exactly by the MoM solver.\n");
			fprintf(stderr, "The system of equations is singular; -X suppresses the singular fallback.\n");
			return 1;
		}
		/* Singular MoM system: solve it exactly through other class
		 * orders, merged stations or classes, or interpolation of the
		 * exact limit of a perturbation (util/singular.h).  This replaces
		 * the earlier automatic perturbation at digit 20, which gave an
		 * approximation and no exact G. */
		bool report = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output;
		int nres = SG_NRES(qnm);
		mpq_t* exact = sg_mpq_arr(nres);
		char tier[256];
		/* Tier 4 runs on mommod's engine: the Ntot+1 perturbed solves of
		 * interpolation done over Q would each be a full exact solve with
		 * long rationals, while per prime in Z/p they are word-size
		 * recursions, combined by CRT into the same exact value. */
		mommod_config(mom_threads, 61, 0, 0, 64);
		sg_ops ops = { sg_direct, mommod_interp, NULL, " per prime" };
		int rc = sg_solve_tiered(qnm, exact, tier, sizeof tier, &ops, 1);
		if (rc != 0) {
			fprintf(stderr, "\nError: the MoM system is singular for this model, and every exact fallback\n"
			                "failed (other class orders, merged stations or classes, interpolation).\n"
			                "The normalizing constant may be zero or undefined; bin/ca solves any\n"
			                "well-posed model by convolution.  -p digit gives a perturbed approximation.\n");
			return 1;
		}
		fprintf(stderr, "mom: singular MoM system, solved exactly via %s\n", tier);
		if (report) printf("\nSingular MoM system, solved exactly via %s\n", tier);
		report_exact(qnm, exact, verbose_output, log_output, normconst_output, normconst_g_output,
		             throughput_output, queue_output, debug_output, scale_factor);
		sg_mpq_arr_free(exact, nres);
	} else if (result == -1) {
		fprintf(stderr, "\nError: Model cannot be solved even with perturbation.\n");
		return 1;
	}
	
	// Clean up original parameter storage
	if(perturbation_digit > 0) {
		for(int i = 0; i < qnm->M; i++) {
			for(int j = 0; j < qnm->R; j++) {
				mpz_clear(original_L[i][j]);
			}
			free(original_L[i]);
		}
		free(original_L);
		for(int j = 0; j < qnm->R; j++) {
			mpz_clear(original_Z[j]);
		}
		free(original_Z);
	}
	
	mpz_clear(scale_factor);
	return 0;
}

