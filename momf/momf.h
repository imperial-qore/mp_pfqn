#ifndef LINEAR
#define LINEAR
//#define MPFR
#include <fpla.h>

typedef struct
{
	int n;
	int k;
	int card; /* nchoosek(n+k-1,k) */
	int **combs; /* matrix of nchoosek(n+k-1,k) combinations with repetition of n elements */
}combsrep;

typedef struct
{
	fp_mat_t diag; /* diagonal block, overwritten in place by its LU factors */
	fp_mat_t orig; /* pristine copy of the block, kept only when refinement is on */
	int      dim;  /* order of the block */
	int*      lu_indices; /* indices returned by fp_ludcmp for diagblock */
	fp_msp_t nondiag; /* non diagonal coefficients */
}fp_diagblock;

typedef struct
{
	int m;
	int r;
	int H;
	int numdiagblocks;
	fp_diagblock* C; 
	fp_msp_t A12;
	fp_msp_t B1r;
	fp_msp_t B2r;
}LS;

#include "profiling.h"
#include "util.h"
extern qnmodel* qnm;

extern bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;

/* Wilkinson/Moler iterative refinement controls.
 * REFINE_ITERS = 0 disables refinement, so that the solver is a plain
 * fixed-precision MoM and the instability of the recursion can be
 * measured in isolation. */
extern int    REFINE_ITERS;   /* max refinement sweeps per block solve */
extern double REFINE_TOL;     /* stopping tolerance on the relative correction */
extern int    REFINE_SWEEPS;  /* running total of sweeps actually performed */
extern int    REFINE_SOLVES;  /* running total of block solves */
extern bool   DIAGNOSE;       /* per-step recursion diagnostics on stderr */

#define copy_Gk_in_g(Gk,g){int t; for(t=0;t<cardGk;t++) fp_set(g[t],Gk[t]); }
#define copy_G_in_g(G,g){int t; for(t=0;t<cardG;t++) fp_set(g[cardGk+t],G[t]); }
#define free_Gk() for (t=0;t<cardGk;t++) fp_clear(Gk[t]); free(Gk);
#define free_G() for (t=0;t<cardG;t++) fp_clear(G[t]); free(G);

fp_vec_t blocksolve(LS* ls, fp_vec_t b);
int mdecrease(qnmodel* qnm, fp_vec_t G, fp_vec_t Gk, fp_vec_t g, fp_vec_t gr, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor);
void perfindices(qnmodel* qnm, fp_vec_t G, fp_vec_t Gk, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor);
LS* setupls(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r);
LS* setupls_progress(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r, double setup_start, int show_progress);
void freels(LS* ls);


#endif
