#ifndef LINEAR
#define LINEAR
//#define MPFR
#include <zpla.h>

typedef struct
{
	int n;
	int k;
	int card; /* nchoosek(n+k-1,k) */
	int **combs; /* matrix of nchoosek(n+k-1,k) combinations with repetition of n elements */
}combsrep;

typedef struct
{
	zp_mat_t diag; /* diagonal block, overwritten in place by its LU factors */
	int      dim;  /* order of the block */
	int*      lu_indices; /* indices returned by zp_ludcmp for diagblock */
	zp_msp_t nondiag; /* non diagonal coefficients */
}zp_diagblock;

typedef struct
{
	int m;
	int r;
	int H;
	int numdiagblocks;
	zp_diagblock* C; 
	zp_msp_t A12;
	zp_msp_t B1r;
	zp_msp_t B2r;
}LS;

#include "profiling.h"
#include "util.h"
extern qnmodel* qnm;

extern bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;

extern int VERBOSE_MOD;   /* progress reporting on stderr */

/* Number of quantities exported from each modular image:
 * G, then X_r for each class, then Q_kr for each station and class. */
#define NRES(qnm) (1 + (qnm)->R + (qnm)->M * (qnm)->R)

#define copy_Gk_in_g(Gk,g){int t; for(t=0;t<cardGk;t++) zp_set(g[t],Gk[t]); }
#define copy_G_in_g(G,g){int t; for(t=0;t<cardG;t++) zp_set(g[cardGk+t],G[t]); }
#define free_Gk() for (t=0;t<cardGk;t++) zp_clear(Gk[t]); free(Gk);
#define free_G() for (t=0;t<cardG;t++) zp_clear(G[t]); free(G);

zp_vec_t blocksolve(LS* ls, zp_vec_t b);
int mdecrease(qnmodel* qnm, zp_vec_t G, zp_vec_t Gk, zp_vec_t g, zp_vec_t gr, zp_t* res);
void exportres(qnmodel* qnm, zp_vec_t G, zp_vec_t Gk, zp_t* res);
LS* setupls(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r);
LS* setupls_progress(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r, double setup_start, int show_progress);
void freels(LS* ls);


#endif
