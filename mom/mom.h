#ifndef LINEAR
#define LINEAR
//#define MPFR
#include <gmpla.h>

typedef struct
{
	int n;
	int k;
	int card; /* nchoosek(n+k-1,k) */
	int **combs; /* matrix of nchoosek(n+k-1,k) combinations with repetition of n elements */
}combsrep;

typedef struct
{
	mpq_mat_t diag; /* diagonal block matrix */
	int*      lu_indices; /* indices returned by mpq_ludcmp for diagblock */
	mpq_msp_t nondiag; /* non diagonal coefficients */
}mpq_diagblock;

typedef struct
{
	int m;
	int r;
	int H;
	int numdiagblocks;
	mpq_diagblock* C; 
	mpq_msp_t A12;
	mpq_msp_t B1r;
	mpq_msp_t B2r;
}LS;

#include "profiling.h"
#include "util.h"
extern qnmodel* qnm;

extern bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;

#define copy_Gk_in_g(Gk,g){int t; for(t=0;t<cardGk;t++) mpq_set(g[t],Gk[t]); }
#define copy_G_in_g(G,g){int t; for(t=0;t<cardG;t++) mpq_set(g[cardGk+t],G[t]); }
#define free_Gk() for (t=0;t<cardGk;t++) mpq_clear(Gk[t]); free(Gk);
#define free_G() for (t=0;t<cardG;t++) mpq_clear(G[t]); free(G);

mpq_vec_t blocksolve(LS* ls, mpq_vec_t b);
int mdecrease(qnmodel* qnm, mpq_vec_t G, mpq_vec_t Gk, mpq_vec_t g, mpq_vec_t gr, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, bool debug_full_precision, mpz_t scale_factor);
void perfindices(qnmodel* qnm, mpq_vec_t G, mpq_vec_t Gk, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor);
LS* setupls(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r);


#endif
