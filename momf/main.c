#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <fpla.h>
#include "momf.h"
#include "lvtab.h"
#include "singular.h"

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
mpq_t* MOM_CAPTURE = NULL;
int MOM_CAPTURE_BAD = 0;

/* Wilkinson/Moler refinement controls (declared in momf.h) */
int    REFINE_ITERS  = 0;
double REFINE_TOL    = 0.0;
int    REFINE_SWEEPS = 0;
int    REFINE_SOLVES = 0;
bool   DIAGNOSE      = false;



/* --- recursion diagnostics ------------------------------------------
 * Reports, per population step, the dynamic range of the unknown vector
 * and the amount of cancellation in forming the right-hand side
 * b = B1r*g - A12*G.  A cancellation ratio of 10^k means k digits are
 * destroyed in that subtraction at every step, which decides whether the
 * instability lives in the linear solve (curable by iterative
 * refinement) or in the recursion itself (not curable). */
static void diagreport(int* n, int R, fp_vec_t g, int ng,
                       fp_vec_t t1v, fp_vec_t t2v, fp_vec_t b, int nb)
{
	int i,s;
	double gmax=0.0, gmin=0.0, amax=0.0, bmax=0.0, v;
	for (i=0;i<ng;i++)
	{
		v = fabs(fp_get_d(g[i]));
		if (v > gmax) gmax = v;
		if (v > 0.0 && (gmin == 0.0 || v < gmin)) gmin = v;
	}
	for (i=0;i<nb;i++)
	{
		v = fabs(fp_get_d(t1v[i])); if (v > amax) amax = v;
		v = fabs(fp_get_d(t2v[i])); if (v > amax) amax = v;
		v = fabs(fp_get_d(b[i]));   if (v > bmax) bmax = v;
	}
	fprintf(stderr, "diag n=(");
	for (s=0;s<R-1;s++) fprintf(stderr,"%d,",n[s]);
	fprintf(stderr, "%d) |g|max=%.3e |g|min=%.3e range=%.3e cancel=%.3e\n",
	        n[R-1], gmax, gmin, (gmin>0.0)?gmax/gmin:0.0,
	        (bmax>0.0)?amax/bmax:0.0);
}

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
static void lv_free(fp_vec_t v, long n)
{
	long t;
	if (!v) return;
	for (t=0;t<=n;t++) fp_clear(v[t]);
	free(v);
}

/* Advance levels 0..r-2 from population n-1_r to n by the class-r
 * population constraint.  glow points at the old level r-1 (second part
 * of g), which must not yet be overwritten. */
static void lower_step(qnmodel* qnm, lvtab* T, int r, fp_vec_t* low, fp_t* glow, long nr)
{
	int l, i, j, s, m = T->m;
	fp_t acc, c, w, zr;
	fp_init(acc); fp_init(c); fp_init(w); fp_init(zr);
	fp_set_z(zr, qnm->Z[r-1]);
	for (l=0;l<=r-2;l++) {
		fp_t* above = (l+1 == r-1) ? glow : low[l+1];
		for (i=0;i<T->card[l];i++) {
			for (s=0;s<r;s++) {
				fp_mul(acc, zr, low[l][i*r+s]);
				for (j=0;j<m;j++) {
					if (mpz_sgn(qnm->L[j][r-1]) == 0) continue;
					fp_set_z(c, qnm->L[j][r-1]);
					fp_set_si(w, qnm->mi[j] + T->cb[l][i][j], 1);
					fp_mul(c, c, w);
					fp_mul(w, c, above[T->up[l][i*m+j]*r+s]);
					fp_add(acc, acc, w);
				}
				fp_set_si(w, nr, 1);
				fp_div(low[l][i*r+s], acc, w);
			}
		}
	}
	fp_clear(acc); fp_clear(c); fp_clear(w); fp_clear(zr);
}

int solve_model(qnmodel* qnm, int argc, char** argv, bool verbose_output, bool log_output, 
                bool normconst_output, bool normconst_g_output, bool throughput_output, 
                bool queue_output, bool debug_output, bool bounds_output, 
                int perturbation_digit, int perturbation_seed, bool auto_perturbation,
                mpz_t scale_factor, mpz_t** original_L, mpz_t* original_Z) {
	int t,i,r,s,Ntot,h,k;
	int *n;
	int cardG, cardGk;
	fp_t tmp; 
	fp_t tmp2; 
	fp_vec_t g=NULL,Gk=NULL,G=NULL,b=NULL,b2=NULL,gr=NULL;
	LS* A=NULL;
	double step_start, step_elapsed;
	double class_start, class_elapsed;
	double setup_start, setup_elapsed;
	
	fp_init(tmp);
	fp_init(tmp2);
	
	/* initializations */
	nckinit(qnm->M+qnm->R-1,qnm->R); /* declare maximum nck term */
	n=(int*)int_vec(qnm->R,0); /* current population vector */
	int m=qnm->M;
	Ntot=0;
	int R=qnm->R, l;
	int full = mom_mixed_think(qnm); /* else the mdecrease descent is exact */
	lvtab* T = full ? lvtab_new(m,R) : NULL;
	fp_vec_t *low  = (fp_vec_t*)calloc(R>1?R-1:1,sizeof(fp_vec_t)); /* levels 0..r-2 */
	fp_vec_t *lowr = (fp_vec_t*)calloc(R>1?R-1:1,sizeof(fp_vec_t)); /* their snapshot at N-1_r */
	
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
			b=(fp_vec_t)fp_vec(m,0,1);
			b2=(fp_vec_t)fp_vec(m,0,1);
			/* the norm.consts for n=(0,0,...,0) are equal to 1 */
			g=(fp_vec_t)fp_vec(m+1,1,1); 
			gr=(fp_vec_t)fp_vec(m+1,0,1); /* G(n-1_1)=0 at n=0: stays 0 if N_1=0 */
			Gk=(fp_vec_t)fp_vec(m,1,1);
			G=(fp_vec_t)fp_vec(1,1,1);
		}
		else
		{
			freels(A);
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
					fp_t* old  = (l == r-2) ? g  + cardGk_old : low[l];
					fp_t* snap = (l == r-2) ? gr + cardGk_old : lowr[l];
					fp_vec_t nw = (fp_vec_t)fp_vec(T->card[l]*r,0,1);
					for (i=0;i<T->card[l];i++) {
						for (s=0;s<=r-2;s++) fp_set(nw[i*r+s], old[i*(r-1)+s]);
						fp_set(nw[i*r+r-1], snap[i*(r-1)]);
					}
					if (l < r-2) { lv_free(low[l],T->card[l]*(r-1)); lv_free(lowr[l],T->card[l]*(r-1)); }
					low[l]  = nw;
					lowr[l] = (fp_vec_t)fp_vec(T->card[l]*r,0,1); /* stays 0 if N_r=0 */
				}
			}
			for (t=0;t<nck(m+r-2,r-1)*(r-1);t++)
			{
				fp_clear(b[t]);
				fp_clear(b2[t]);
			}
			for (t=0;t<nck(m+r-3,r-2)*(r-1);t++)
				fp_clear(G[t]);
			b=(fp_vec_t)fp_vec(cardGk,0,1);
			b2=(fp_vec_t)fp_vec(cardGk,0,1);
			G=(fp_vec_t)fp_vec(cardG,0,1);
			Gk=(fp_vec_t)fp_vec(cardGk,0,1);
			for (i=0;i<nck(m+r-2,r-1);i++)
			{
				for (s=0;s<=r-2;s++)
				{
					fp_set(G[i*r+s],g[i*(r-1)+s]);
				}
				fp_set(G[i*r+r-1],gr[i*(r-1)]);
			}
			fp_mspvecmul(b,A->A12,G);
			for(t=0;t<cardGk;t++)
				fp_neg(b[t],b[t]);	
			Gk=(fp_vec_t)blocksolve(A,b);
			if (Gk == NULL) {
				return -1; // Signal singular matrix
			}
			for (t=0;t<nck(m+r-2,r-1)*(r-1);t++)
			{
				fp_clear(g[t]);
				fp_clear(gr[t]);
			}
			g=(fp_vec_t)fp_vec(cardGk+cardG,0,1);
			gr=(fp_vec_t)fp_vec(cardGk+cardG,0,1);
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
				fp_vecdup(gr,g,cardGk+cardG);
				for (l=0;full && l<=r-2;l++) fp_vecdup(lowr[l],low[l],T->card[l]*r);
			}

			/* compute G=B2r*gr/nr */

			fp_mspvecmul(G,A->B2r,g);
			fp_t ninv; fp_init(ninv); fp_set_si(ninv,1,n[r-1]); 
			for(t=0;t<cardG;t++) fp_mul(G[t],ninv,G[t]);
			if (full) lower_step(qnm,T,r,low,g+cardGk,n[r-1]);

			/* compute b=B1r*gr - A12*G */
			fp_mspvecmul(b,A->B1r,g);
			fp_mspvecmul(b2,A->A12,G);
			if (DIAGNOSE)
			{
				fp_vec_t bkeep=(fp_vec_t)fp_vec(cardGk,0,1);
				fp_vecdup(bkeep,b,cardGk);
				for(t=0;t<cardGk;t++) fp_sub(b[t],b[t],b2[t]);
				diagreport(n,qnm->R,g,cardGk+cardG,bkeep,b2,b,cardGk);
				for(t=0;t<=cardGk;t++) fp_clear(bkeep[t]);
				free(bkeep);
			}
			else
			for(t=0;t<cardGk;t++) fp_sub(b[t],b[t],b2[t]);

			/* compute Gk */
			Gk=(fp_vec_t)blocksolve(A,b);
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
		fp_t *L0, *S0, *L1, *S1;
		int cardGkR = nck(m+R-1,R)*R;
		L0 = (R >= 2) ? low[0]  : g + cardGkR;
		S0 = (R >= 2) ? lowr[0] : gr + cardGkR;
		L1 = (R >= 3) ? low[1]  : (R == 2 ? g + cardGkR  : g);
		S1 = (R >= 3) ? lowr[1] : (R == 2 ? gr + cardGkR : gr);
		fp_vec_t G0 = (fp_vec_t)fp_vec(R+1,0,1);
		fp_vec_t G1 = (fp_vec_t)fp_vec(m*(R+1),0,1);
		for (s=0;s<R;s++) fp_set(G0[s], L0[s]);
		fp_set(G0[R], S0[0]);
		for (k=0;k<m;k++) {
			for (s=0;s<R;s++) fp_set(G1[k*(R+1)+s], L1[k*R+s]);
			fp_set(G1[k*(R+1)+R], S1[k*R]);
		}
		perfindices(qnm,G0,G1,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
		lv_free(G0,R+1); lv_free(G1,m*(R+1));
	}
	if (full) for (l=0;l<R-1;l++) { lv_free(low[l],T->card[l]*R); lv_free(lowr[l],T->card[l]*R); }
	free(low); free(lowr);
	lvtab_free(T);
	fp_clear(tmp);
	fp_clear(tmp2);
	
	return 0; // Success
}

/* Solve of a derived model for the tiered fallback (util/singular.h): the
 * constants are captured in out[], converted exactly from the working
 * precision, instead of printed.  0 on success, 1 if the MoM system of q
 * is singular or the result is not finite. */
static int sg_direct(qnmodel* q, mpq_t* out, void* ctx)
{
	mpz_t one; int rc;
	(void)ctx;
	mpz_init_set_ui(one, 1);
	MOM_CAPTURE = out; MOM_CAPTURE_BAD = 0;
	rc = solve_model(q, 0, NULL, false, false, false, false, false, true, false, false,
	                 0, 0, false, one, NULL, NULL);
	MOM_CAPTURE = NULL;
	mpz_clear(one);
	return (rc == 0 && !MOM_CAPTURE_BAD) ? 0 : 1;
}

/* A singular system is singular in floating point only when a pivot
 * rounds to exactly zero, which rounding usually prevents: the solve then
 * "succeeds" with garbage scaled by the reciprocal of a rounding error.
 * Tiers 1-3 therefore reject relative LU pivots below 2^-(P/2), and accept
 * a solve only if a second solve at about twice the working precision
 * (same pivot threshold) agrees with it to half the working precision;
 * otherwise it counts as singular and the next tier is tried.  Entries
 * that vanish identically (N_r = 0) are set to zero exactly. */
static int sg_direct_checked(qnmodel* q, mpq_t* out, void* ctx)
{
	int nres = SG_NRES(q), u, k, rc;
	mpfr_prec_t w0 = FP_WPREC, r0 = FP_RPREC;
	mpq_t* hi = sg_mpq_arr(nres);
	int piv0 = FP_PIVTOL_BITS;
	if (FP_PIVTOL_BITS < (int)(w0/2)) FP_PIVTOL_BITS = (int)(w0/2);
	rc = sg_direct(q, out, ctx);
	if (rc == 0)
	{
		FP_WPREC = 2*w0 + 64;
		FP_RPREC = 2*r0 + 64;
		if (FP_RPREC < FP_WPREC) FP_RPREC = FP_WPREC;
		mpfr_set_default_prec(FP_WPREC);
		rc = sg_direct(q, hi, ctx);
		FP_WPREC = w0; FP_RPREC = r0;
		mpfr_set_default_prec(w0);
	}
	FP_PIVTOL_BITS = piv0;
	for (u=0; rc==0 && u<q->R; u++) if (q->N[u] == 0)
	{
		mpq_set_ui(out[1+u], 0, 1); mpq_set_ui(hi[1+u], 0, 1);
		for (k=0;k<q->M;k++) { mpq_set_ui(out[SG_IDX_GK(q->R,k,u)], 0, 1); mpq_set_ui(hi[SG_IDX_GK(q->R,k,u)], 0, 1); }
	}
	if (rc == 0)
	{
		mpq_t d, lim; mpq_init(d); mpq_init(lim);
		for (u=0; u<nres; u++)
		{
			mpq_sub(d, out[u], hi[u]); mpq_abs(d, d);
			mpq_abs(lim, hi[u]);
			mpq_div_2exp(lim, lim, (mp_bitcnt_t)(w0/2));
			if (mpq_cmp(d, lim) > 0) { rc = 1; break; }
		}
		mpq_clear(d); mpq_clear(lim);
	}
	sg_mpq_arr_free(hi, nres);
	return rc;
}

/* Tier 4 in floating point.  Two errors reach the extrapolated value:
 * the cancellation of the Lagrange sum at t=0, which amplifies the error
 * of the perturbed solves by 2^amp (amp from sg_interp_exact_amp), and
 * the fixed-precision drift of the MoM recursion itself.  The perturbed
 * models are therefore solved at the working precision P plus extra bits,
 * and the result is accepted only once a run at more extra bits agrees
 * with the previous one to P bits and amp is covered by the extra, so
 * the value returned is accurate to about the working precision. */
static int SG_EXTRA_BITS = 0;
static int sg_interp_round(qnmodel* q, mpq_t* out, long extra, double* amp)
{
	mpfr_prec_t w0 = FP_WPREC, r0 = FP_RPREC;
	sg_ops inner = { sg_direct, NULL, NULL, NULL };
	int rc, piv0 = FP_PIVTOL_BITS;
	if (piv0 > 0) FP_PIVTOL_BITS = (int)((w0 + extra)/2);   /* the threshold follows the precision */
	FP_WPREC = w0 + extra;
	FP_RPREC = r0 + extra;
	if (FP_RPREC < FP_WPREC) FP_RPREC = FP_WPREC;
	mpfr_set_default_prec(FP_WPREC);
	rc = sg_interp_exact_amp(q, out, &inner, amp);
	FP_PIVTOL_BITS = piv0;
	FP_WPREC = w0; FP_RPREC = r0;
	mpfr_set_default_prec(w0);
	return rc;
}

static int sg_agree(const mpq_t* a, const mpq_t* b, int n, long bits)
{
	int u, ok = 1;
	mpq_t d, lim; mpq_init(d); mpq_init(lim);
	for (u=0; u<n && ok; u++)
	{
		mpq_sub(d, a[u], b[u]); mpq_abs(d, d);
		mpq_abs(lim, b[u]); mpq_div_2exp(lim, lim, (mp_bitcnt_t)bits);
		if (mpq_cmp(d, lim) > 0) ok = 0;
	}
	mpq_clear(d); mpq_clear(lim);
	return ok;
}

static int sg_interp_fp(qnmodel* q, mpq_t* out, void* ctx)
{
	long P = (long)FP_WPREC, extra = 64;
	int nres = SG_NRES(q), it, rc;
	double amp = 0.0, amp2 = 0.0;
	mpq_t* prev = sg_mpq_arr(nres);
	(void)ctx;
	rc = sg_interp_round(q, prev, extra, &amp);
	for (it = 0; rc == 0 && it < 8; it++)
	{
		long next = 2*extra + P;
		if (next < (long)amp + 64) next = (long)amp + 64;
		rc = sg_interp_round(q, out, next, &amp2);
		if (rc != 0) break;
		if (amp2 + 16 <= (double)next && sg_agree(prev, out, nres, P)) { extra = next; break; }
		extra = next; amp = amp2;
		{ int u; for (u=0;u<nres;u++) mpq_set(prev[u], out[u]); }
		if (it == 7) rc = 1;
	}
	SG_EXTRA_BITS = (int)extra;
	sg_mpq_arr_free(prev, nres);
	return rc;
}

/* Print the reconstructed constants through perfindices, in the layout it
 * expects (see mom/main.c). */
static void report_exact(qnmodel* q, mpq_t* ex, bool verbose_output, bool log_output,
                         bool normconst_output, bool normconst_g_output, bool throughput_output,
                         bool queue_output, bool debug_output, mpz_t scale_factor)
{
	int R = q->R, M = q->M, i, k;
	fp_vec_t G0 = (fp_vec_t)fp_vec(R+1,0,1);
	fp_vec_t G1 = (fp_vec_t)fp_vec(M*(R+1),0,1);
	for (i=0;i<=R;i++) mpfr_set_q(G0[i], ex[i], MPFR_RNDN);
	for (k=0;k<M;k++) for (i=0;i<R;i++) mpfr_set_q(G1[k*(R+1)+1+i], ex[SG_IDX_GK(R,k,i)], MPFR_RNDN);
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
	long working_prec = 53;    /* IEEE double mantissa by default */
	long residual_prec = 0;    /* 0 -> twice the working precision */
	int  refine_iters = 0;     /* 0 -> plain fixed-precision MoM */
	double refine_tol = 0.0;   /* 0 -> always run the full sweep budget */
	long pivtol_bits = -1;     /* -1 -> half the working precision, 0 -> off */

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
		printf("  -w, --prec bits  : Working floating-point precision (default: 53, i.e. IEEE double)\n");
		printf("  --resprec bits   : Residual precision for iterative refinement (default: 2*prec)\n");
		printf("  -I n             : Wilkinson/Moler refinement sweeps per block solve (default: 0 = off)\n");
		printf("  --pivtol bits    : Treat LU pivots below 2^-bits (relative to their row) as singular\n");
		printf("                     (default: P/2; 0 = exact zero only, the earlier behavior)\n");
		printf("  -T tol           : Stop refining when the relative correction falls below tol (default: 0)\n");
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
			printf("  -w, --prec bits  : Working floating-point precision (default: 53, i.e. IEEE double)\n");
			printf("  --resprec bits   : Residual precision for iterative refinement (default: 2*prec)\n");
			printf("  -I n             : Wilkinson/Moler refinement sweeps per block solve (default: 0 = off)\n");
		printf("  --pivtol bits    : Treat LU pivots below 2^-bits (relative to their row) as singular\n");
		printf("                     (default: P/2; 0 = exact zero only, the earlier behavior)\n");
			printf("  -T tol           : Stop refining when the relative correction falls below tol (default: 0)\n");
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
		} else if(strcmp(argv[i], "--pivtol") == 0) {
			if(i + 1 < argc) {
				pivtol_bits = atol(argv[++i]);
			} else {
				printf("Error: --pivtol option requires a bit count\n");
				return -1;
			}
		} else if(strcmp(argv[i], "-D") == 0 || strcmp(argv[i], "--diag") == 0) {
			DIAGNOSE = true;
		} else if(strcmp(argv[i], "-P") == 0 || strcmp(argv[i], "-W") == 0) {
			/* -P means probabilities in promom/procomom/camarg, -W the
			 * witness primes in mommod/comommod */
			fprintf(stderr, "momf: %s is now %s\n", argv[i],
			        argv[i][1] == 'P' ? "-w/--prec" : "--resprec");
			return 1;
		} else if(strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--prec") == 0) {
			if(i + 1 < argc) {
				working_prec = atol(argv[++i]);
				if(working_prec < 2) {
					printf("Error: working precision must be at least 2 bits\n");
					return -1;
				}
			} else {
				printf("Error: -w option requires a bit count\n");
				return -1;
			}
		} else if(strcmp(argv[i], "--resprec") == 0) {
			if(i + 1 < argc) {
				residual_prec = atol(argv[++i]);
			} else {
				printf("Error: --resprec option requires a bit count\n");
				return -1;
			}
		} else if(strcmp(argv[i], "-I") == 0 || strcmp(argv[i], "--refine") == 0) {
			if(i + 1 < argc) {
				refine_iters = atoi(argv[++i]);
			} else {
				printf("Error: -I option requires an iteration count\n");
				return -1;
			}
		} else if(strcmp(argv[i], "-T") == 0 || strcmp(argv[i], "--reftol") == 0) {
			if(i + 1 < argc) {
				refine_tol = atof(argv[++i]);
			} else {
				printf("Error: -T option requires a tolerance\n");
				return -1;
			}
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
		printf("  -w, --prec bits  : Working floating-point precision (default: 53, i.e. IEEE double)\n");
		printf("  --resprec bits   : Residual precision for iterative refinement (default: 2*prec)\n");
		printf("  -I n             : Wilkinson/Moler refinement sweeps per block solve (default: 0 = off)\n");
		printf("  --pivtol bits    : Treat LU pivots below 2^-bits (relative to their row) as singular\n");
		printf("                     (default: P/2; 0 = exact zero only, the earlier behavior)\n");
		printf("  -T tol           : Stop refining when the relative correction falls below tol (default: 0)\n");
		return -1;
	}
	
	// Validate bounds option
	if(bounds_output && perturbation_digit == 0) {
		printf("Error: -b/--bounds option requires -p option to be specified\n");
		return -1;
	}
	
	/* Install the arithmetic precision before anything allocates an
	 * mpfr_t: fp_init sizes every value from FP_WPREC. */
	FP_WPREC = (mpfr_prec_t) working_prec;
	/* A singular system shows up in floating point as a pivot made of
	 * rounding error rather than an exact zero; reject relative pivots
	 * below 2^-(P/2) so that such a model reaches the singular-model
	 * fallback instead of returning garbage. */
	FP_PIVTOL_BITS = (int)(pivtol_bits < 0 ? working_prec/2 : pivtol_bits);
	FP_RPREC = (mpfr_prec_t) (residual_prec > 0 ? residual_prec : 2*working_prec);
	if (FP_RPREC < FP_WPREC) FP_RPREC = FP_WPREC;
	mpfr_set_default_prec(FP_WPREC);
	REFINE_ITERS = refine_iters;
	REFINE_TOL   = refine_tol;

	qnm=(qnmodel*)readmodel(model_file);
	qn_require_closed(qnm, "momf", 0);
	
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
	
	// Try to solve the model
	int result = solve_model(qnm, argc, argv, verbose_output, log_output,
	                        normconst_output, normconst_g_output, throughput_output,
	                        queue_output, debug_output, bounds_output,
	                        perturbation_digit, perturbation_seed, auto_perturbation,
	                        scale_factor, original_L, original_Z);

	// If singular and no perturbation yet, apply automatic perturbation
	if (result == -1 && !auto_perturbation && perturbation_digit == 0) {
		/* Singular MoM system: other class orders, merged stations or
		 * classes, or interpolation of the limit of a perturbation
		 * (util/singular.h), in place of the earlier automatic
		 * perturbation at digit 20.  Tiers 1-3 are exact transformations;
		 * tier 4 runs at raised precision (sg_interp_fp). */
		bool report = !log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output;
		int nres = SG_NRES(qnm);
		mpq_t* exact = sg_mpq_arr(nres);
		char tier[256];
		sg_ops ops = { sg_direct_checked, sg_interp_fp, NULL, NULL };
		int rc = sg_solve_tiered(qnm, exact, tier, sizeof tier, &ops, 1);
		if (rc != 0) {
			fprintf(stderr, "\nError: the MoM system is singular for this model, and every fallback\n"
			                "failed (other class orders, merged stations or classes, interpolation).\n"
			                "The normalizing constant may be zero or undefined; bin/ca solves any\n"
			                "well-posed model by convolution.  -p digit gives a perturbed approximation.\n");
			return 1;
		}
		if (SG_EXTRA_BITS && strncmp(tier, "interpolation", 13) == 0)
			fprintf(stderr, "momf: singular MoM system, solved via %s at %ld+%d bits\n", tier, (long)FP_WPREC, SG_EXTRA_BITS);
		else
			fprintf(stderr, "momf: singular MoM system, solved via %s\n", tier);
		if (report) printf("\nSingular MoM system, solved via %s\n", tier);
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

