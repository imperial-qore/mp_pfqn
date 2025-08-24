#ifndef LINEAR
#define LINEAR
#define MPFR
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
	mpq_mat_t A11;
	mpq_msp_t A12;
	mpq_msp_t B1;
	mpq_msp_t B2;
}LS;

#include "profiling.h"
#include "util.h"
extern qnmodel* qnm;

#define copy_Gk_in_g(Gk,g){int t; for(t=0;t<cardGk;t++) mpq_set(g[t],Gk[t]); }
#define copy_G_in_g(G,g){int t; for(t=0;t<cardG;t++) mpq_set(g[cardGk+t],G[t]); }
#define free_Gk() for (t=0;t<cardGk;t++) mpq_clear(Gk[t]); free(Gk);
#define free_G() for (t=0;t<cardG;t++) mpq_clear(G[t]); free(G);

LS* setupls(combsrep *Dn, qnmodel* qnm, int* N, int r);
int hash(combsrep* Dn, int *comb, int i);
#endif
