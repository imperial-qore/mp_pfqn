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

/* BTF block decomposition types */
typedef struct {
	int size;           /* dimension of square diagonal block */
	mpq_mat_t diag;     /* dense sub-matrix [size x size], holds LU factors after factorize */
	int* lu_indices;    /* LU pivot indices, NULL before factorization */
	int level;          /* number of non-zeros in pattern (popcount of bitmask) */
} btf_block;

typedef struct {
	int from_block;     /* row block index (lower level, source of CE/PC row) */
	int to_block;       /* column block index (higher level, destination of coupling) */
	int nrows, ncols;
	mpq_mat_t mat;      /* dense sub-matrix [nrows x ncols] */
} btf_offdiag;

typedef struct {
	int num_blocks;
	btf_block* blocks;
	int* block_starts;      /* block_starts[b] = first permuted index of block b */
	int num_offdiag;
	btf_offdiag* offdiag;
	int* row_perm;          /* row_perm[permuted_i] = original A11 row */
	int* col_perm;          /* col_perm[permuted_i] = original A11 col */
	int total_size;         /* = cardGk */
} btf_info;

typedef struct
{
	int m;
	int r;
	mpq_mat_t A11;
	mpq_msp_t A12;
	mpq_msp_t B1;
	mpq_msp_t B2;
	btf_info* btf;
	/* Direct LU factorization (fallback when BTF has singular blocks) */
	mpq_mat_t A11_lu;       /* LU-factored copy of A11 (NULL if BTF works) */
	int* A11_lu_indices;    /* LU pivot indices for direct A11 solve */
	/* Reduced system for non-filtered combo handling */
	int n_filt_cols;            /* number of filtered columns (0 = not used) */
	int* filt_col_map;          /* filt_col_map[j] = original A11 col for filtered col j */
	mpq_mat_t A11_reduced;      /* LU-factored square matrix [n_filt_cols x n_filt_cols] */
	int* A11_reduced_indices;   /* LU pivot indices from mpq_ludcmp */
	int* pivot_rows;            /* original A11 row indices of independent rows [n_filt_cols] */
	/* Non-filtered combo info (for direct Gk computation at solve time) */
	int n_nf_combos;            /* number of non-filtered combos */
	int* nf_combo_indices;      /* combo indices in Dn that are non-filtered [n_nf_combos] */
	int* nf_col_indices;        /* A11 col indices for non-filtered Gk [n_nf_combos * M] */
}LS;

#include "profiling.h"
#include "util.h"
extern qnmodel* qnm;

#define copy_Gk_in_g(Gk,g){int t; for(t=0;t<cardGk;t++) mpq_set(g[t],Gk[t]); }
#define copy_G_in_g(G,g){int t; for(t=0;t<cardG;t++) mpq_set(g[cardGk+t],G[t]); }
#define free_Gk() for (t=0;t<cardGk;t++) mpq_clear(Gk[t]); free(Gk);
#define free_G() for (t=0;t<cardG;t++) mpq_clear(G[t]); free(G);

LS* setupls(combsrep *Dn, qnmodel* qnm, int* N, int r, double setup_start, int show_progress);
int hash(combsrep* Dn, int *comb, int i);
int* mpq_ludcmp_progress(mpq_mat_t A, int N, double setup_start, int* n_vec, int R, int show_progress);

/* BTF block decomposition functions */
btf_info* btf_decompose(mpq_mat_t A11, combsrep* Dn, int r, int M, int cardGk);
int btf_factorize(btf_info* btf);
int btf_solve(btf_info* btf, mpq_vec_t b, int cardGk);
void btf_free(btf_info* btf);
#endif
