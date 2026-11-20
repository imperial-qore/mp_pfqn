#ifndef COMOMMOD_H
#define COMOMMOD_H

/* comommod - multi-modular CoMoM.
 *
 * The exact rational CoMoM solver (bin/comom) carries normalizing
 * constants whose length grows as Theta(N log N) bits, so every
 * arithmetic operation in the class-oriented recursion gets more
 * expensive as the recursion proceeds.  Here the whole recursion is
 * replayed independently in Z/p for a sequence of word-size primes and
 * the rational answers are recovered at the end by incremental CRT plus
 * rational reconstruction, exactly as mommod does for MoM.
 *
 * Every source in this directory is a transliteration of its comom/
 * counterpart with mpq_t replaced by a residue in a machine word; the
 * row order, column order and coefficients are unchanged, so the two
 * solvers walk the same recursion and must agree.
 */

#include <zpla.h>

typedef struct
{
	int n;          /* length of a combination */
	int k;          /* = M */
	int card;       /* nchoosek(n+k-1,k) */
	int **combs;    /* the combinations, sorted by sortbynnzpos */
	/* Optional lookup accelerating hash().  comb -> row by binary search
	 * on a mixed-radix encoding, built once per class and shared by every
	 * modular image; NULL when the encoding would overflow, in which case
	 * hash() falls back to the linear scan comom uses. */
	unsigned long long *keys;   /* keys[i] = code of combs[order[i]] */
	int *order;
} combsrep;

/* BTF block decomposition, over Z/p */
typedef struct {
	int size;
	zp_mat_t diag;      /* holds the LU factors after btf_factorize */
	int* lu_indices;
	int level;
} btf_block;

typedef struct {
	int from_block;
	int to_block;
	int nrows, ncols;
	zp_mat_t mat;
} btf_offdiag;

typedef struct {
	int num_blocks;
	btf_block* blocks;
	int* block_starts;
	int num_offdiag;
	btf_offdiag* offdiag;
	int* row_perm;
	int* col_perm;
	int total_size;
} btf_info;

typedef struct
{
	int m;
	int r;
	long int cardG;
	long int cardGk;
	zp_mat_t A11;
	zp_msp_t A12;
	zp_msp_t B1;
	zp_msp_t B2;
	btf_info* btf;
	zp_mat_t A11_lu;        /* LU of A11, used when BTF is singular mod p */
	int* A11_lu_indices;
}LS;

#include "profiling.h"
#include "util.h"

extern qnmodel* qnm;
extern bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;
extern int VERBOSE_MOD;

/* Quantities exported from each modular image:
 *   res[0]                        = G(N)
 *   res[r]                        = G(N-e_r)          (numerator of X_r)
 *   res[1 + R + (k-1)*R + (r-1)]  = G^k(N-e_r)        (numerator of Q_kr)
 * The raw constants are exported rather than the ratios: X_r and Q_kr are
 * rationals whose numerator and denominator are each about as long as G,
 * so reconstructing them directly would need roughly four times the
 * modulus, hence four times as many primes.  The divisions are done once,
 * exactly, after reconstruction. */
#define NRES(qnm) (1 + (qnm)->R + (qnm)->M * (qnm)->R)

combsrep* makeDn(int r, int M);
void      freeDn(combsrep* Dn);
int       hash(combsrep* Dn, int *comb, int i);

LS*  setupls(combsrep *Dn, qnmodel* qnm, int *N, int r);
void freels(LS* ls);

btf_info* btf_decompose(zp_mat_t A11, combsrep* Dn, int r, int M, int cardGk);
int  btf_factorize(btf_info* btf);
int  btf_solve(btf_info* btf, zp_vec_t b, int cardGk);
void btf_free(btf_info* btf);

int  comom_image(qnmodel* qnm, combsrep** Dn_all, zp_t p, zp_t* res);

#endif
