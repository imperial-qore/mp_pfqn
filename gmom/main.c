#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"
#include "lvtab.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static int NTHREADS = 1;   /* -j: sub-models with the same queue count in parallel */
/* multiplicities of the model being solved: station k has mi[k]/mden copies
 * (util.h, qn_mult), split in the model tree into a never-removed base
 * (mi[k] % mden)/mden plus mi[k]/mden removable unit copies */
static const int* GM_MI = NULL;
static int GM_MDEN = 1;

struct rusage ruse;
double t0, t1;

extern void gmva(mpq_t G, mpz_t** dem, int nq, int* pop, mpz_t* Z, int r);
extern void gmom_descend(qnmodel* qn, mpq_vec_t vlM, mpq_vec_t grM, int ltop, mpq_vec_t* Gout, mpq_vec_t* Gkout);
extern void gmom_measures(qnmodel* qn, mpq_vec_t G, mpq_vec_t Gk, int oe, int og, int ol, int ot, int oq, int* cperm, mpz_t scale_factor, gmom_out* cap);

/* ---- small dense helpers over rectangular mpq matrices ---------------- */

/* rop[nrows] = A[nrows x ncols] * v[ncols] */
static void matvec_rect(mpq_vec_t rop, mpq_mat_t A, int nrows, int ncols, mpq_vec_t v)
{
	int i, j; mpq_t t; mpq_init(t);
	for (i = 0; i < nrows; i++) {
		mpq_set_si(rop[i], 0, 1);
		for (j = 0; j < ncols; j++)
			if (mpq_sgn(A[i][j]) && mpq_sgn(v[j])) { mpq_mul(t, A[i][j], v[j]); mpq_add(rop[i], rop[i], t); }
	}
	mpq_clear(t);
}

/* rop[ncols] = A^T[ncols x nrows] * v[nrows] */
static void tvec_rect(mpq_vec_t rop, mpq_mat_t A, int nrows, int ncols, mpq_vec_t v)
{
	int i, j; mpq_t t; mpq_init(t);
	for (j = 0; j < ncols; j++) {
		mpq_set_si(rop[j], 0, 1);
		for (i = 0; i < nrows; i++)
			if (mpq_sgn(A[i][j]) && mpq_sgn(v[i])) { mpq_mul(t, A[i][j], v[i]); mpq_add(rop[j], rop[j], t); }
	}
	mpq_clear(t);
}

/* LU factorization of A^T A, cached across a class sweep.  A (setup1's
 * coefficient matrix) holds demands only: the populations and think times
 * sit in B, C, D, E and F.  So A^T A is the same at every step of the
 * class-r sweep and is factored once per (level, class), as Sec. 5 of the
 * paper assumes, instead of once per population step. */
typedef struct { mpq_mat_t AtA; int* indx; int ncols; } nfact;

static void nfact_free(nfact* f);

/* may_fail: return NULL on a singular A^T A instead of exiting */
static nfact* nfact_new(mpq_mat_t A, int nrows, int ncols, int may_fail)
{
	nfact* f = (nfact*) malloc(sizeof(nfact));
	f->ncols = ncols;
	f->AtA = mpq_matzeros(ncols, ncols);
	mpq_mattransmul(f->AtA, A, A, nrows, ncols);
	f->indx = mpq_ludcmp(f->AtA, ncols);
	if (!f->indx && may_fail) { f->indx = (int*) malloc(sizeof(int)); nfact_free(f); return NULL; }
	if (!f->indx) { fprintf(stderr, "gmom: SINGULAR A^T A  nrows=%d ncols=%d\n", nrows, ncols);
		if (getenv("GMOM_DBG")) { fprintf(stderr,"A =\n"); mpq_matprint(A, nrows, ncols); }
		exit(2); }
	return f;
}

static void nfact_free(nfact* f)
{
	int i, j;
	if (!f) return;
	for (i = 0; i < f->ncols+1; i++) { for (j = 0; j < f->ncols+1; j++) mpq_clear(f->AtA[i][j]); free(f->AtA[i]); }
	free(f->AtA); free(f->indx); free(f);
}

/* x = (A^T A)^{-1} A^T rhs with the cached factorization; fresh x[ncols] */
static mpq_vec_t nfact_solve(nfact* f, mpq_mat_t A, int nrows, mpq_vec_t rhs)
{
	mpq_vec_t x = mpq_vec(f->ncols, 0, 1);
	tvec_rect(x, A, nrows, f->ncols, rhs);
	mpq_lubksb(f->AtA, x, f->ncols, f->indx);
	return x;
}

/* G of the sub-model on stations q[0..d-1], station q[j] with multiplicity
 * en[j]/GM_MDEN plus combo[j] replicas, at population pop, think times Z,
 * r classes.  Integer multiplicities go to gmva as repeated demand rows. */
static void base_Gq(mpq_t G, mpz_t** L, const int* q, int d, const int* en, int* combo, int r, int* pop, mpz_t* Z)
{
	int j, c, nq = 0;
	if (GM_MDEN > 1) {
		mpz_t** dem = (mpz_t**) malloc(d * sizeof(mpz_t*));
		int* e = (int*) malloc(d * sizeof(int));
		for (j = 0; j < d; j++) { dem[j] = L[q[j]]; e[j] = en[j] + GM_MDEN * combo[j]; }
		gbase_exp(G, dem, e, GM_MDEN, d, pop, Z, r);
		free(dem); free(e);
		return;
	}
	for (j = 0; j < d; j++) nq += en[j] + combo[j];
	mpz_t** dem = (mpz_t**) malloc(nq * sizeof(mpz_t*));
	c = 0;
	for (j = 0; j < d; j++) { int t; for (t = 0; t < en[j] + combo[j]; t++) dem[c++] = L[q[j]]; }
	gmva(G, dem, nq, pop, Z, r);
	free(dem);
}

/* the full model (all m stations, multiplicities GM_MI) with combo replicas */
static void base_G(mpq_t G, mpz_t** L, int m, int* combo, int r, int* pop, mpz_t* Z)
{
	int* q = (int*) malloc(m * sizeof(int)), j;
	for (j = 0; j < m; j++) q[j] = j;
	base_Gq(G, L, q, m, GM_MI, combo, r, pop, Z);
	free(q);
}

static void oner(int* pop, int* nvec, int r, int s)
{
	int c; for (c = 0; c < r; c++) pop[c] = nvec[c];
	if (s >= 1) pop[s-1] -= 1;
}

/* free/alloc an mpq vector slot */
static mpq_vec_t vnew(int n) { return mpq_vec(n, 0, 1); }
static void vfree(mpq_vec_t v, int n) { int i; if (!v) return; for (i = 0; i < n+1; i++) mpq_clear(v[i]); free(v); }
static void fp_copy_vec(mpq_vec_t d, mpq_vec_t s, int n) { int i; for (i = 0; i < n; i++) mpq_set(d[i], s[i]); }

/* Lower replica levels 0..ltop-2 of the full M-queue basis (see
 * util/lvtab.h).  The full model's two slots hold levels ltop (vl) and
 * ltop-1 (vl1), ltop = r-1 under --bmax 1; lv1 = ltop-1 below.  The
 * levels below are carried for the full model only, filled exactly by
 * convolution when class r starts and advanced along the class-r sweep
 * by the class-r population constraint
 *   n_r G^{+d}(n) = Z_r G^{+d}(n-1_r) + sum_j L_jr (m_j+d_j) G^{+d+1_j}(n-1_r),
 * which holds for any think time.  Levels 0 and 1 at N then give G(N),
 * X and Q without a downward pass. */
static void fill_low(lvtab* T, mpq_vec_t* low, mpz_t** L, int r, int* nvec, mpz_t* Z, int lv1)
{
	int l, i, s, M = T->m;
	int* pop = (int*) malloc(r * sizeof(int));
	for (l = 0; l <= lv1-1; l++)
		for (i = 0; i < T->card[l]; i++)
			for (s = 0; s < r; s++) {
				oner(pop, nvec, r, s);
				base_G(low[l][i*r+s], L, M, T->cb[l][i], r, pop, Z);
			}
	free(pop);
}

/* advance levels 0..r-3 from n-1_r to n; above2 is the old level r-2 */
static void lower_step(lvtab* T, int* mi, mpz_t** L, mpz_t* Z, int r, mpq_vec_t* low, mpq_vec_t above2, int nr, int lv1)
{
	int l, i, j, s, M = T->m;
	mpq_t acc, c, w, zr;
	mpq_init(acc); mpq_init(c); mpq_init(w); mpq_init(zr);
	mpq_set_z(zr, Z[r-1]);
	for (l = 0; l <= lv1-1; l++) {
		mpq_vec_t above = (l+1 == lv1) ? above2 : low[l+1];
		for (i = 0; i < T->card[l]; i++)
			for (s = 0; s < r; s++) {
				mpq_mul(acc, zr, low[l][i*r+s]);
				for (j = 0; j < M; j++) {
					if (mpz_sgn(L[j][r-1]) == 0) continue;
					mpq_set_z(c, L[j][r-1]);
					mpq_set_si(w, mi[j] + (long)GM_MDEN * T->cb[l][i][j], GM_MDEN);
					mpq_canonicalize(w);
					mpq_mul(c, c, w);
					mpq_mul(w, c, above[T->up[l][i*M+j]*r+s]);
					mpq_add(acc, acc, w);
				}
				mpq_set_si(w, nr, 1);
				mpq_div(low[l][i*r+s], acc, w);
			}
	}
	mpq_clear(acc); mpq_clear(c); mpq_clear(w); mpq_clear(zr);
}


/* ---- model tree of the variable-branching recursion (--bmax) --------- */

/* A sub-model is a copy vector v: station k has the base (mi[k] % mden)/mden,
 * never removed, plus v_k unit copies, so its multiplicity is en/mden with
 * en = mi % mden + mden v_k (util.h, qn_mult); its support q are the
 * stations with a positive multiplicity.  It has ntag = min{Bmax, #removable}
 * GCEs, a removable station being one with a copy whose removal leaves a
 * non-empty model, and basis level l = max{1, r - Bc}, Bc = ntag (Bc = 1 for a
 * single-station support, which keeps the original single-queue level)
 * while r classes are processed (Theorem 2).  Removing a copy changes the
 * number of removable stations by at most one, so every child carries
 * level l or l+1, the two slots a GCE row reads.  With unit multiplicities
 * this is the original tree: v is the queue subset, removable = every queue
 * of a model with d >= 2, Bc = min{d, Bmax}. */
typedef struct {
	int d, ntag, Bc, tot;
	int* v;          /* copy vector over all M stations */
	int* q;          /* support: original station indices, increasing */
	int* en;         /* multiplicity numerators of q[], over GM_MDEN */
	int* tagpos;     /* child c removes a copy of local station tagpos[c] */
	int* cdrop;      /* ... and drops it from its support */
	int* child;
	int* clevel;     /* child basis levels, current class */
	int l;           /* basis level, current class */
	mpq_vec_t vl, vl1; int lenl, lenl1;
	mpq_vec_t pvl, pvl1; int lenp, lenp1;   /* both at n-1_r */
	LSB* ls; nfact* fa;
} gmodel;

typedef struct { gmodel* m; int n; } gtree;

static int gm_lev(const gmodel* g, int r) { int v = r - g->Bc; return v < 1 ? 1 : v; }

/* m[0] is the full model (v_k = floor(mi_k)); each model's children remove
 * one copy of each of its last ntag removable stations; children shared
 * between parents are built once.  A child has one copy fewer, so models
 * are appended by decreasing total copies tot and the array read backwards
 * is a bottom-up order.  Unit multiplicities, Bmax = 1: the prefix chain
 * {1..M} -> {1..M-1} -> ... -> {1}. */
static gtree* gtree_new(int M, int bmax, const int* mi, int mden)
{
	gtree* T = (gtree*) calloc(1, sizeof(gtree));
	int cap = 16, i, c, a, j, k;
	T->m = (gmodel*) calloc(cap, sizeof(gmodel));
	T->m[0].v = (int*) malloc(M * sizeof(int));
	for (k = 0; k < M; k++) T->m[0].v[k] = mi[k] / mden;
	T->n = 1;
	int* cv = (int*) malloc(M * sizeof(int));
	for (i = 0; i < T->n; i++) {
		gmodel* g = &T->m[i];
		int d = 0, tot = 0, nrem = 0;
		for (k = 0; k < M; k++) { if (g->v[k] > 0 || mi[k] % mden) d++; tot += g->v[k]; }
		g->d = d; g->tot = tot;
		g->q = (int*) malloc((d ? d : 1) * sizeof(int)); g->en = (int*) malloc((d ? d : 1) * sizeof(int));
		for (k = 0, a = 0; k < M; k++) if (g->v[k] > 0 || mi[k] % mden) { g->q[a] = k; g->en[a] = mi[k] % mden + mden * g->v[k]; a++; }
		int* rem = (int*) malloc((d ? d : 1) * sizeof(int));
		for (a = 0; a < d; a++) { k = g->q[a]; if (g->v[k] >= 1 && (d >= 2 || g->v[k] >= 2 || mi[k] % mden)) rem[nrem++] = a; }
		int ntag = nrem < bmax ? nrem : bmax;
		g->ntag = ntag; g->Bc = (d == 1) ? 1 : ntag;
		g->child  = (int*) calloc(ntag > 0 ? ntag : 1, sizeof(int));
		g->clevel = (int*) calloc(ntag > 0 ? ntag : 1, sizeof(int));
		g->tagpos = (int*) calloc(ntag > 0 ? ntag : 1, sizeof(int));
		g->cdrop  = (int*) calloc(ntag > 0 ? ntag : 1, sizeof(int));
		for (c = 0; c < ntag; c++) {
			a = rem[nrem - ntag + c]; k = T->m[i].q[a];
			memcpy(cv, T->m[i].v, M * sizeof(int)); cv[k]--;
			T->m[i].tagpos[c] = a;
			T->m[i].cdrop[c] = (cv[k] == 0 && mi[k] % mden == 0);
			int found = -1;
			for (j = i+1; j < T->n && found < 0; j++)
				if (T->m[j].tot == tot-1 && !memcmp(T->m[j].v, cv, M * sizeof(int))) found = j;
			if (found < 0) {
				if (T->n == cap) { cap *= 2; T->m = (gmodel*) realloc(T->m, cap * sizeof(gmodel)); }
				memset(&T->m[T->n], 0, sizeof(gmodel));
				T->m[T->n].v = (int*) malloc(M * sizeof(int)); memcpy(T->m[T->n].v, cv, M * sizeof(int));
				T->m[T->n].tot = tot-1;
				found = T->n++;
			}
			T->m[i].child[c] = found;
		}
		free(rem);
	}
	free(cv);
	return T;
}

static void gtree_free(gtree* T)
{
	int i;
	if (!T) return;
	for (i = 0; i < T->n; i++) {
		gmodel* g = &T->m[i];
		free(g->q); free(g->child); free(g->clevel); free(g->v); free(g->en); free(g->tagpos); free(g->cdrop);
		vfree(g->vl, g->lenl); vfree(g->vl1, g->lenl1);
		vfree(g->pvl, g->lenp); vfree(g->pvl1, g->lenp1);
		lsb_free(g->ls); nfact_free(g->fa);
	}
	free(T->m); free(T);
}

/* r-component copy of an (r-1)-component basis: entry i keeps its r-1
 * values and gains G(n-1_{r-1}) = component 0 of the snapshot pv */
static mpq_vec_t vext(mpq_vec_t v, mpq_vec_t pv, int ncomb, int r, int* len)
{
	int i, s;
	mpq_vec_t w = vnew(ncomb * r);
	for (i = 0; i < ncomb; i++) {
		for (s = 0; s < r-1; s++) mpq_set(w[i*r+s], v[i*(r-1)+s]);
		mpq_set(w[i*r+r-1], pv[i*(r-1)]);
	}
	*len = ncomb * r;
	return w;
}

/* new level-l basis of g from rhs = B*vprev + D*vl1 + child GCE terms;
 * vprev = NULL is the class start, where G(n-1_r) = 0 */
static mpq_vec_t level_solve(gmodel* g, gtree* GT, mpq_vec_t vprev)
{
	LSB* ls = g->ls;
	int s, nrows = ls->nrows;
	mpq_vec_t rd = vnew(nrows), rhs = vnew(nrows);
	if (vprev) matvec_rect(rhs, ls->B, nrows, ls->colA, vprev);
	matvec_rect(rd, ls->D, nrows, ls->colD, g->vl1);
	for (s = 0; s < nrows; s++) {
		mpq_add(rhs[s], rhs[s], rd[s]);
		if (ls->cchild[s] >= 0) {
			gmodel* ch = &GT->m[g->child[ls->cchild[s]]];
			mpq_add(rhs[s], rhs[s], (ls->cslot[s] ? ch->vl1 : ch->vl)[ls->ccol[s]]);
		}
	}
	mpq_vec_t x = nfact_solve(g->fa, ls->A, nrows, rhs);
	vfree(rd, nrows); vfree(rhs, nrows);
	return x;
}

/* fill a level-lv basis of model g by exact convolution at nvec */
static mpq_vec_t conv_level(gmodel* g, mpz_t** L, int lv, int r, int* nvec, mpz_t* Z, int* len)
{
	int i, s, sz = nck(g->d + lv - 1, lv);
	int** I = sortbynnzpos(multichoose(g->d, lv), sz, g->d);
	int* pop = (int*) malloc(r * sizeof(int));
	mpq_vec_t v = vnew(sz * r);
	for (i = 0; i < sz; i++) for (s = 0; s < r; s++) { oner(pop, nvec, r, s); base_Gq(v[i*r+s], L, g->q, g->d, g->en, I[i], r, pop, Z); }
	for (i = 0; i < sz; i++) free(I[i]);
	free(I); free(pop);
	*len = sz * r;
	return v;
}

/* both slots of g (levels g->l and g->l-1) by exact convolution: the
 * class initialisation, and every step for the single-queue models */
static void fill_model(gmodel* g, mpz_t** L, int r, int* nvec, mpz_t* Z)
{
	vfree(g->vl, g->lenl);   g->vl  = conv_level(g, L, g->l,   r, nvec, Z, &g->lenl);
	vfree(g->vl1, g->lenl1); g->vl1 = conv_level(g, L, g->l-1, r, nvec, Z, &g->lenl1);
}

/* Singularity screen for --bmax > 1 (pfqn_recursion_singular builds the
 * Bmax = 1 matrices only): every level matrix is loading-only, so test
 * A^T A of every sub-model for every class count once at start. */
static int gb_singular(mpz_t** L, int M, int R, int bmax)
{
	int r, i, c, sing = 0;
	gtree* T = gtree_new(M, bmax, GM_MI, GM_MDEN);
	for (r = 2; r <= R && !sing; r++)
		for (i = 0; i < T->n && !sing; i++) {
			gmodel* g = &T->m[i];
			if (g->d < 2) continue;
			int l = gm_lev(g, r);
			for (c = 0; c < g->ntag; c++) g->clevel[c] = gm_lev(&T->m[g->child[c]], r);
			LSB* ls = setupb(L, g->q, g->d, g->en, GM_MDEN, g->ntag, g->tagpos, g->cdrop, g->clevel, r, l, NULL, NULL, 1);
			mpq_mat_t AtA = mpq_matzeros(ls->colA, ls->colA);
			mpq_mattransmul(AtA, ls->A, ls->A, ls->nrows, ls->colA);
			int* idx = mpq_ludcmp(AtA, ls->colA);
			if (!idx) sing = 1; else free(idx);
			int a, b;
			for (a = 0; a < ls->colA+1; a++) { for (b = 0; b < ls->colA+1; b++) mpq_clear(AtA[a][b]); free(AtA[a]); }
			free(AtA);
			lsb_free(ls);
		}
	gtree_free(T);
	return sing;
}

/* class to move to the recursion position so that the screen passes, or
 * -1 (same contract as pfqn_nonsingular_recclass) */
static int gb_nonsingular_recclass(mpz_t** L, int M, int R, int bmax)
{
	int c, k, s, j;
	if (!gb_singular(L, M, R, bmax)) return R-1;
	for (c = 0; c < R-1; c++) {
		mpz_t** P = (mpz_t**) malloc(M * sizeof(mpz_t*));
		for (k = 0; k < M; k++) {
			P[k] = (mpz_t*) malloc(R * sizeof(mpz_t));
			for (s = 0, j = 0; s < R; s++) if (s != c) mpz_init_set(P[k][j++], L[k][s]);
			mpz_init_set(P[k][R-1], L[k][c]);
		}
		int sing = gb_singular(P, M, R, bmax);
		for (k = 0; k < M; k++) { for (s = 0; s < R; s++) mpz_clear(P[k][s]); free(P[k]); }
		free(P);
		if (!sing) return c;
	}
	return -1;
}

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

/* sign = +1 adds the perturbation (mom's convention), sign = -1 subtracts it.
 * A demand of 0 is left at 0 under sign = -1, since a negative demand is not
 * a model; that entry then carries no perturbation error either way. */
void apply_signed_perturbation_to_model(qnmodel* qnm, int perturbation_digit, int perturbation_seed, int sign, mpz_t scale_factor) {
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
			int zero = (mpz_sgn(qnm->L[i][j]) == 0);
			mpz_mul(qnm->L[i][j], qnm->L[i][j], scale_factor);
			if (sign > 0)          mpz_add_ui(qnm->L[i][j], qnm->L[i][j], perm[i]);
			else if (!zero)        mpz_sub_ui(qnm->L[i][j], qnm->L[i][j], perm[i]);
		}

		// Apply perturbation to Z[j] using the last element of permutation
		int zeroz = (mpz_sgn(qnm->Z[j]) == 0);
		mpz_mul(qnm->Z[j], qnm->Z[j], scale_factor);
		if (sign > 0)             mpz_add_ui(qnm->Z[j], qnm->Z[j], perm[qnm->M]);
		else if (!zeroz)          mpz_sub_ui(qnm->Z[j], qnm->Z[j], perm[qnm->M]);
	}

	free(perm);
}

void apply_perturbation_to_model(qnmodel* qnm, int perturbation_digit, int perturbation_seed, mpz_t scale_factor) {
	apply_signed_perturbation_to_model(qnm, perturbation_digit, perturbation_seed, +1, scale_factor);
}

/* A completed gmom solve: everything downstream of the recursion (measures,
 * validation) works off this, so -b can run the recursion twice under two
 * perturbation seeds and bracket the two answers. */
typedef struct {
	qnmodel*  qn;
	mpq_vec_t G0;     /* level 0 at N: G(N), G(N-1_s)            */
	mpq_vec_t G1;     /* level 1 at N: G(M+1_k, N), G(M+1_k, N-1_s) */
	int*      cperm;
	mpz_t     scale;
} gmom_sol;

/* One complete solve: model load, domain checks, singularity screen and the
 * divide-and-conquer recursion.  Returns 0 on success, 1 on error, 2 when
 * --validate consumed the run. */
static int gmom_compute(char* model_file, int pdigit, int pseed, int psign,
                        bool out_e, bool force_perturb, bool no_perturb, bool validate,
                        bool print_model, int bmax, gmom_sol* sol)
{
	int i, s, r, nr;
	qn_allow_rational_mi = 1;   /* fractional multiplicities (util.h, qn_mult) */
	qnmodel* qn = readmodel(model_file);
	qn_require_closed(qn, "gmom", 0);
	int M = qn->M, R = qn->R;
	if (R < 2) { fprintf(stderr, "gmom requires at least two classes\n"); return 1; }
	{ int c; for (c = 0; c < R; c++) if (qn->N[c] < 0) {
		fprintf(stderr, "gmom: open/mixed models are not supported (class %d is open); use a closed model\n", c+1); return 1; } }
	if (qn->hasOpen) { fprintf(stderr, "gmom: open/mixed models are not supported\n"); return 1; }
	if (qn->isLD)    { fprintf(stderr, "gmom: load-dependent stations are not supported\n"); return 1; }
	/* multiplicities (integer copies or rational, util.h qn_mult) live in
	 * the model tree; the shared Bmax = 1 oracle assumes unit ones */
	GM_MI = qn->mi; GM_MDEN = qn_mden(qn);
	int unit_mi = 1;
	{ int q; for (q = 0; q < M; q++) if (qn->mi[q] != GM_MDEN) unit_mi = 0; }
	mpz_t** L = qn->L; mpz_t* Z = qn->Z; int* N = qn->N;
	nckinit(M + R - 1, R);
	/* nck() caches into a global table on first use, which is not thread
	 * safe: fill every entry the recursion reads before going parallel */
	{ int a, b; for (a = 0; a <= M + R + 1 && a < MAXNCKTABLE; a++) for (b = 0; b <= a; b++) (void)nck(a, b); }

	/* keep the unperturbed demands so the header can report the model as
	 * read plus the per-class epsilon actually applied (as mom does) */
	mpz_t** original_L = (mpz_t**) malloc(M * sizeof(mpz_t*));
	mpz_t*  original_Z = (mpz_t*)  malloc(R * sizeof(mpz_t));
	{ int k, j;
	  for (k = 0; k < M; k++) { original_L[k] = (mpz_t*) malloc(R * sizeof(mpz_t));
	    for (j = 0; j < R; j++) mpz_init_set(original_L[k][j], qn->L[k][j]); }
	  for (j = 0; j < R; j++) mpz_init_set(original_Z[j], qn->Z[j]); }
	int* original_N = (int*) malloc(R * sizeof(int));
	{ int j; for (j = 0; j < R; j++) original_N[j] = qn->N[j]; }

	mpz_t scale_factor; mpz_init_set_ui(scale_factor, 1);
	int applied_digit = 0;
	if (pdigit > 0) { apply_signed_perturbation_to_model(qn, pdigit, pseed, psign, scale_factor); applied_digit = pdigit; }

	/* Init-time singularity screen (loading-only, one factorisation per
	 * level).  If the recursion would hit a degenerate matrix, move the
	 * offending class to the recursion position (its loadings then leave
	 * the coefficient matrix); if no single-class reorder removes it, the
	 * model is genuinely degenerate for gmom's basis. */
	int* cperm = (int*) malloc(R * sizeof(int));
	for (i = 0; i < R; i++) cperm[i] = i;
	/* Bmax = 1 uses the shared oracle (same matrices); larger Bmax adds
	 * GCE rows and changes the levels, so it screens its own matrices */
	if (bmax == 1 && unit_mi ? pfqn_recursion_singular(L, M, R) : gb_singular(L, M, R, bmax)) {
		int fix = bmax == 1 && unit_mi ? pfqn_nonsingular_recclass(L, M, R) : gb_nonsingular_recclass(L, M, R, bmax);
		if (fix < 0) {
			if (no_perturb) {
				/* exact-only probe (safe_gmom): report and fail rather
				 * than degrade to a perturbed approximation */
				fprintf(stderr, "gmom: singular system under every class order and --no-perturb given; no exact answer\n");
				return 1;
			}
			if (out_e && !force_perturb) {
				fprintf(stderr, "gmom: model is singular under every class order (genuine loading degeneracy);\n  no exact G available. Use bin/ca or bin/safe_gmom, or -g/-t/-q with -p for a perturbed\n  approximation, or --force-perturb to print the perturbed rational anyway.\n");
				return 1;
			}
			/* auto-perturb at digit 20 (as mom/comom do) for an approximate answer */
			if (mpz_cmp_ui(scale_factor, 1) == 0) {
				fprintf(stderr, "gmom: singular system; automatically applying perturbation at digit 20 (approximate result).\n");
				apply_perturbation_to_model(qn, 20, pseed, scale_factor);
				applied_digit = 20;
				L = qn->L; Z = qn->Z; N = qn->N;
			}
			if (out_e) fprintf(stderr, "gmom: --force-perturb: the rational printed under -e is the EXACT G of the\n  PERTURBED model, not of the model as read.\n");
		} else {
		/* cperm[j] = original class placed at position j (fix goes last) */
		int c = 0, j;
		for (j = 0; j < R; j++) if (j != fix) cperm[c++] = j;
		cperm[R-1] = fix;
		/* apply the permutation to L columns, N and Z (work in place) */
		mpz_t** L2 = (mpz_t**) malloc(M * sizeof(mpz_t*));
		int k;
		for (k = 0; k < M; k++) { L2[k] = (mpz_t*) malloc(R * sizeof(mpz_t)); for (j = 0; j < R; j++) { mpz_init_set(L2[k][j], L[k][cperm[j]]); } }
		int* N2 = (int*) malloc(R * sizeof(int)); mpz_t* Z2 = (mpz_t*) malloc(R * sizeof(mpz_t));
		for (j = 0; j < R; j++) { N2[j] = N[cperm[j]]; mpz_init_set(Z2[j], Z[cperm[j]]); }
		L = L2; N = N2; Z = Z2;
		qn->L = L2; qn->N = N2; qn->Z = Z2;   /* measures reads from qn */
		}
	}

	/* Header in the default (no output flag) mode.  Printed off a view of
	 * the model in its ORIGINAL class order and with the demands as read,
	 * so a perturbed run is self-describing: the eps= column reports the
	 * perturbation actually applied, whether from -p or from the automatic
	 * digit-20 fallback above. */
	if (print_model) {
		qnmodel view = *qn;
		view.L = original_L; view.Z = original_Z; view.N = original_N;
		if (applied_digit > 0)
			printmodel_with_perturbation(&view, applied_digit, scale_factor, pseed, original_L, original_Z);
		else
			printmodel(&view);
	}

	int* nvec = (int*) calloc(R, sizeof(int));

	gtree* GT = gtree_new(M, bmax, GM_MI, GM_MDEN);
	gmodel* top = &GT->m[0];
	mpq_vec_t grM = NULL; int lengrM = 0;   /* top vl at N-1_R */
	mpq_vec_t gr1M = NULL; int lengr1M = 0; /* top vl1 at N-1_R */
	int full = mom_mixed_think(qn); /* else the descent (gmom_descend) is exact */
	lvtab* T = full ? lvtab_new(M, R) : NULL;
	mpq_vec_t* low  = (mpq_vec_t*) calloc(R, sizeof(mpq_vec_t)); /* levels 0..ltop-2, top model */
	mpq_vec_t* lowr = (mpq_vec_t*) calloc(R, sizeof(mpq_vec_t)); /* same at N-1_R */
	int* lowlen = (int*) calloc(R, sizeof(int));
	int ltop = 1;

	nvec[0] = N[0];
	for (r = 2; r <= R; r++) {
		int ltop_old = ltop;
		/* New class, class r empty.  Levels and level systems, factored
		 * once for the whole sweep.  Class 2 starts by convolution (only
		 * class 1 is populated, so the lattice is one-dimensional).  From
		 * class 3 on, each basis is carried over exactly: the r-1 entries
		 * per combination are the old ones, and the new decrement
		 * G(n-1_{r-1}) is the snapshot taken one step before the end of
		 * the class-(r-1) sweep.  When the level rises (l(d,r) = l(d,r-1)+1)
		 * the old level becomes vl1 and the new vl solves the level system
		 * with the class-r term G(n-1_r) = 0.  Exact for any Z (the
		 * reference's shift initialisation is Z=0-only). */
		if (full) {
			/* lower levels of the full model: old levels carried over, plus
			 * the old vl1 when the full model's level rises */
			int ll, ln = gm_lev(top, r);
			mpq_vec_t* nlow = (mpq_vec_t*) calloc(R, sizeof(mpq_vec_t));
			int* nlen = (int*) calloc(R, sizeof(int));
			if (r == 2) {
				for (ll = 0; ll <= ln-2; ll++) { nlen[ll] = T->card[ll]*r; nlow[ll] = vnew(nlen[ll]); }
				fill_low(T, nlow, L, r, nvec, Z, ln-1);
			} else {
				for (ll = 0; ll <= ln-2; ll++) {
					if (ll <= ltop_old-2) nlow[ll] = vext(low[ll], lowr[ll], T->card[ll], r, &nlen[ll]);
					else                  nlow[ll] = vext(top->vl1, top->pvl1, T->card[ll], r, &nlen[ll]);
				}
			}
			for (ll = 0; ll < R; ll++) { vfree(low[ll], lowlen[ll]); vfree(lowr[ll], lowlen[ll]); lowr[ll] = NULL; }
			free(low); free(lowlen);
			low = nlow; lowlen = nlen;
			for (ll = 0; ll <= ln-2; ll++) lowr[ll] = vnew(lowlen[ll]);
		}
		/* setup and factorization read no other model: all in parallel */
		int* lold = (int*) calloc(GT->n, sizeof(int));
		#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
		for (i = GT->n - 1; i >= 0; i--) {
			gmodel* g = &GT->m[i];
			int c;
			lold[i] = g->l;
			g->l = gm_lev(g, r);
			for (c = 0; c < g->ntag; c++) {
				g->clevel[c] = gm_lev(&GT->m[g->child[c]], r);
				if (g->clevel[c] != g->l && g->clevel[c] != g->l + 1) { fprintf(stderr, "gmom: internal error, child level %d under level %d\n", g->clevel[c], g->l); exit(1); }
			}
			lsb_free(g->ls); g->ls = NULL;
			nfact_free(g->fa); g->fa = NULL;
			/* single-queue models too: without GCEs their system has r rows
			 * and r unknowns at any level, square, and it replaces a lattice
			 * convolution per entry and step (Sec. 4.2).  It is singular when
			 * the queue has a zero demand in a class s < r; those models stay
			 * on convolution (ls = NULL). */
			g->ls = setupb(L, g->q, g->d, g->en, GM_MDEN, g->ntag, g->tagpos, g->cdrop, g->clevel, r, g->l, nvec, Z, 0);
			g->fa = nfact_new(g->ls->A, g->ls->nrows, g->ls->colA, g->d == 1);
			if (!g->fa) { lsb_free(g->ls); g->ls = NULL; }
		}
		/* class start: a level solve reads the children (d-1 queues), so
		 * one group of equal d at a time, bottom-up, each group in parallel */
		for (int gb = GT->n; gb > 0; ) {
		int ga = gb - 1;
		while (ga > 0 && GT->m[ga-1].tot == GT->m[gb-1].tot) ga--;
		#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
		for (i = gb - 1; i >= ga; i--) {
			gmodel* g = &GT->m[i];
			int lo = lold[i];
			if (r == 2 || !g->ls) {
				fill_model(g, L, r, nvec, Z);
			} else if (g->l == lo) {
				int n1, n2;
				mpq_vec_t a = vext(g->vl,  g->pvl,  g->lenl/(r-1),  r, &n1);
				mpq_vec_t b = vext(g->vl1, g->pvl1, g->lenl1/(r-1), r, &n2);
				vfree(g->vl, g->lenl); vfree(g->vl1, g->lenl1);
				g->vl = a; g->lenl = n1; g->vl1 = b; g->lenl1 = n2;
			} else {
				int n1;
				mpq_vec_t b = vext(g->vl, g->pvl, g->lenl/(r-1), r, &n1);
				vfree(g->vl1, g->lenl1); g->vl1 = b; g->lenl1 = n1;
				vfree(g->vl, g->lenl);
				g->vl = level_solve(g, GT, NULL); g->lenl = g->ls->colA;
			}
		}
		gb = ga;
		}
		free(lold);
		ltop = top->l;
		/* G(n-1_r) snapshots: zero until the last step of the sweep sets
		 * them (they stay zero when N_r = 0) */
		for (i = 0; i < GT->n; i++) {
			gmodel* g = &GT->m[i];
			vfree(g->pvl, g->lenp);  g->pvl  = vnew(g->lenl);  g->lenp  = g->lenl;
			vfree(g->pvl1, g->lenp1); g->pvl1 = vnew(g->lenl1); g->lenp1 = g->lenl1;
		}

		/* population sweep of class r */
		for (nr = 1; nr <= N[r-1]; nr++) {
			nvec[r-1] = nr;
			if (nr == N[r-1]) {   /* snapshot every basis at n-1_r */
				int ll;
				for (i = 0; i < GT->n; i++) {
					gmodel* g = &GT->m[i];
					fp_copy_vec(g->pvl, g->vl, g->lenl);
					fp_copy_vec(g->pvl1, g->vl1, g->lenl1);
				}
				for (ll = 0; full && ll <= ltop-2; ll++) fp_copy_vec(lowr[ll], low[ll], lowlen[ll]);
			}
			if (full) lower_step(T, qn->mi, L, Z, r, low, top->vl1, nr, ltop-1);
			for (int gb = GT->n; gb > 0; ) {   /* bottom-up by queue count, each group in parallel */
			int ga = gb - 1;
			while (ga > 0 && GT->m[ga-1].tot == GT->m[gb-1].tot) ga--;
			#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
			for (i = gb - 1; i >= ga; i--) {
				gmodel* g = &GT->m[i];
				if (!g->ls) { fill_model(g, L, r, nvec, Z); continue; }
				LSB* ls = g->ls;
				int colA = ls->colA, colD = ls->colD;
				/* level l-1: (E*vl + F*vl1) / nr, from the previous step */
				mpq_vec_t ue = vnew(colD), uf = vnew(colD), u = vnew(colD);
				matvec_rect(ue, ls->E, colD, colA, g->vl);
				matvec_rect(uf, ls->F, colD, colD, g->vl1);
				mpq_t inv; mpq_init(inv); mpq_set_si(inv, 1, nr);
				for (int s = 0; s < colD; s++) { mpq_add(u[s], ue[s], uf[s]); mpq_mul(u[s], u[s], inv); }
				mpq_clear(inv);
				vfree(g->vl1, g->lenl1); g->vl1 = u; g->lenl1 = colD;
				vfree(ue, colD); vfree(uf, colD);
				mpq_vec_t x = level_solve(g, GT, g->vl);
				vfree(g->vl, g->lenl); g->vl = x; g->lenl = colA;
			}
			gb = ga;
			}
		}
	}
	/* the full model at N-1_R */
	grM = top->pvl; lengrM = top->lenp; top->pvl = NULL;
	gr1M = top->pvl1; lengr1M = top->lenp1; top->pvl1 = NULL;

	/* ===================== output / validation ===================== */
	if (validate) {
		mpq_t ref; mpq_init(ref);
		int total = 0, totbad = 0, sl;
		for (i = 0; i < GT->n; i++) {
			gmodel* g = &GT->m[i];
			int badm = 0;
			for (sl = 0; sl < 2; sl++) {
				int len, k;
				mpq_vec_t v = sl ? g->vl1 : g->vl;
				mpq_vec_t c = conv_level(g, L, sl ? g->l-1 : g->l, R, N, Z, &len);
				for (k = 0; k < len; k++) {
					total++;
					if (mpq_cmp(c[k], v[k]) != 0) {
						badm++; totbad++;
						if (badm <= 2) {
							int a;
							fprintf(stderr, "  model {");
							for (a = 0; a < g->d; a++) fprintf(stderr, "%s%d", a ? "," : "", g->q[a]+1);
							gmp_fprintf(stderr, "} level %d entry %d: gmom=%Qd conv=%Qd\n", sl ? g->l-1 : g->l, k, v[k], c[k]);
						}
					}
				}
				vfree(c, len);
			}
		}
		mpq_clear(ref);
		if (totbad) { printf("VALIDATE: %d/%d WRONG\n", totbad, total); return 1; }
		printf("VALIDATE: all %d basis entries of %d sub-models match exact convolution\n", total, GT->n);
		return 2;
	}

	/* G(N), X, Q come from levels 0 and 1 of the full basis at N.
	 * N_R=0: the class-R sweep never runs, and G(.,N-1_R)=0 (negative
	 * population), so the snapshots stay zero. */
	if (!grM || !gr1M) { fprintf(stderr, "gmom: internal error, N-1_R basis not captured\n"); return 1; }
	if (!full) {
		gmom_descend(qn, top->vl, grM, ltop, &sol->G0, &sol->G1);
	} else {
		/* level lv lives in top vl (ltop), top vl1 (ltop-1) or low[lv] */
		mpq_vec_t L0 = (ltop == 1) ? top->vl1 : low[0];
		mpq_vec_t S0 = (ltop == 1) ? gr1M     : lowr[0];
		mpq_vec_t L1 = (ltop == 1) ? top->vl  : (ltop == 2 ? top->vl1 : low[1]);
		mpq_vec_t S1 = (ltop == 1) ? grM      : (ltop == 2 ? gr1M     : lowr[1]);
		mpq_vec_t G0 = vnew(R+1), G1 = vnew(M*(R+1));
		int k;
		for (s = 0; s < R; s++) mpq_set(G0[s], L0[s]);
		mpq_set(G0[R], S0[0]);
		for (k = 0; k < M; k++) {
			for (s = 0; s < R; s++) mpq_set(G1[k*(R+1)+s], L1[k*R+s]);
			mpq_set(G1[k*(R+1)+R], S1[k*R]);
		}
		sol->G0 = G0; sol->G1 = G1;
	}
	{
		int ll;
		for (ll = 0; ll < R; ll++) { vfree(low[ll], lowlen[ll]); vfree(lowr[ll], lowlen[ll]); }
		free(low); free(lowr); free(lowlen);
		vfree(grM, lengrM);
		vfree(gr1M, lengr1M);
		lvtab_free(T);
		gtree_free(GT);
	}
	sol->qn    = qn;
	sol->cperm = cperm;
	mpz_init_set(sol->scale, scale_factor);
	return 0;
}

/* -b: two solves at the same digit and seed, one with +eps and one with -eps.
 * G is a polynomial with non-negative coefficients in the demands, so it is
 * monotone in them and [G(-eps), G(+eps)] encloses the exact G(N) (likewise
 * log G).  X and Q are ratios of such polynomials and are NOT monotone in
 * general, so their interval is an error indicator, not a proven enclosure;
 * the output labels which is which. */
static void print_bracket(const char* label, double a, double b)
{
	double lo = a < b ? a : b, hi = a < b ? b : a;
	double mid = 0.5*(lo+hi);
	double rel = (mid != 0.0) ? (hi-lo)/(2.0*(mid < 0 ? -mid : mid)) : 0.0;
	printf("%-12s %.15e  [%.15e, %.15e]  relhw=%.2e\n", label, mid, lo, hi, rel);
}

int main(int argc, char** argv)
{
	int i;
	bool out_e = false, out_g = false, out_l = false, out_t = false, out_q = false;
	bool validate = false, out_b = false, force_perturb = false, no_perturb = false;
	int pdigit = 0, pseed = 23000, bmax = 1;
	char* model_file = NULL;

	t0 = CPUTIME;
	for (i = 1; i < argc; i++) {
		if      (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--ex"))  out_e = true;
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--nc"))  out_g = true;
		else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--log")) out_l = true;
		else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tput")) out_t = true;
		else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) out_q = true;
		else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bounds")) out_b = true;
		else if (!strcmp(argv[i], "--force-perturb")) force_perturb = true;
		else if (!strcmp(argv[i], "-X") || !strcmp(argv[i], "--no-perturb")) no_perturb = true;
		else if (!strcmp(argv[i], "--validate")) validate = true;
		else if (!strcmp(argv[i], "--bmax")) { if (i+1 < argc) bmax = atoi(argv[++i]); }
		else if (!strcmp(argv[i], "-j")) { if (i+1 < argc) NTHREADS = atoi(argv[++i]); if (NTHREADS < 1) NTHREADS = 1; }
		else if (!strcmp(argv[i], "-p")) { if (i+1 < argc) pdigit = atoi(argv[++i]); }
		else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--seed")) { if (i+1 < argc) pseed = atoi(argv[++i]); }
		else if (argv[i][0] != '-') model_file = argv[i];
	}
	if (!model_file) {
		printf("USAGE: %s [-e|-g|-l|-t|-q] [-b] [--bmax B] [-j n] [--validate] [--force-perturb] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  gmom: generalized (divide-and-conquer) Method of Moments, MoM(Bmax).\n");
		printf("  Outputs the top normalizing-constant basis entry V{M,l}(1),\n");
		printf("  matching the reference mbmom1.  --validate checks the whole\n");
		printf("  basis against exact convolution of the augmented models.\n");
		printf("  --bmax B          : maximum branching factor (default 1): a model with d\n");
		printf("                      queues adds min(d,B) convolution equations and carries\n");
		printf("                      basis level max(1,R-min(d,B)); B >= R-1 gives level 1\n");
		printf("  -j n              : threads, models with the same queue count in parallel\n");
		printf("  -b, --bounds      : bracket a perturbed run by re-solving with a\n");
		printf("                      second seed; requires -p (see gmom/README.md)\n");
		printf("  --force-perturb   : allow the perturbed solve under -e (the printed\n");
		printf("                      rational is then exact for the PERTURBED model)\n");
		printf("  -X, --no-perturb  : fail on a singular system instead of auto-perturbing\n");
		printf("                      (exact-only probe, used by bin/safe_gmom)\n");
		return -1;
	}
	if (bmax < 1) { fprintf(stderr, "gmom: --bmax must be >= 1 (Bmax = 0 is the original MoM: use bin/mom)\n"); return 1; }
	if (out_b && pdigit == 0) {
		fprintf(stderr, "gmom: -b/--bounds requires -p digit (there is nothing to bracket without a perturbation)\n");
		return 1;
	}
	if (out_b && out_e) {
		fprintf(stderr, "gmom: -b/--bounds is incompatible with -e (a bracket is not a single exact rational)\n");
		return 1;
	}

	if (!out_b) {
		gmom_sol sol;
		bool print_model = !(out_e || out_g || out_l || out_t || out_q || validate);
		int rc = gmom_compute(model_file, pdigit, pseed, +1, out_e, force_perturb, no_perturb, validate, print_model, bmax, &sol);
		if (rc) return rc == 2 ? 0 : rc;
		gmom_measures(sol.qn, sol.G0, sol.G1, out_e, out_g, out_l, out_t, out_q, sol.cperm, sol.scale, NULL);
		return 0;
	}

	/* bounds: same digit and seed, one solve with +eps and one with -eps */
	gmom_sol sp, sm;
	int rc = gmom_compute(model_file, pdigit, pseed, +1, false, force_perturb, no_perturb, false, false, bmax, &sp);
	if (rc) return rc == 2 ? 0 : rc;
	int M = sp.qn->M, R = sp.qn->R;
	gmom_out op, om;
	op.X = (double*) malloc(R * sizeof(double)); op.Q = (double*) malloc(M*R * sizeof(double));
	om.X = (double*) malloc(R * sizeof(double)); om.Q = (double*) malloc(M*R * sizeof(double));
	gmom_measures(sp.qn, sp.G0, sp.G1, 0, 0, 0, 0, 0, sp.cperm, sp.scale, &op);
	rc = gmom_compute(model_file, pdigit, pseed, -1, false, force_perturb, no_perturb, false, false, bmax, &sm);
	if (rc) { fprintf(stderr, "gmom: the -eps solve failed; no bracket available\n"); return rc == 2 ? 0 : rc; }
	gmom_measures(sm.qn, sm.G0, sm.G1, 0, 0, 0, 0, 0, sm.cperm, sm.scale, &om);

	printf("========== gmom (perturbation bracket, digit %d, seed %d) ==========\n", pdigit, pseed);
	printf("Columns: midpoint  [lower, upper]  relative half-width.\n");
	printf("Lower/upper are the -eps and +eps solves.  G and log(G) are monotone\n");
	printf("in the demands, so their interval encloses the exact value; the X and Q\n");
	printf("intervals are ratios of such quantities and are error indicators only.\n\n");
	if (out_l) {
		print_bracket("log(G)", om.logG, op.logG);
	} else if (out_g || !(out_t || out_q)) {
		print_bracket("G", om.G, op.G);
	}
	if (out_t || !(out_g || out_l || out_q)) {
		int r;
		printf("\nX (throughputs, indicative interval):\n");
		for (r = 0; r < R; r++) { char lbl[32]; snprintf(lbl, sizeof lbl, "X[%d]", r+1); print_bracket(lbl, om.X[r], op.X[r]); }
	}
	if (out_q || !(out_g || out_l || out_t)) {
		int k, r;
		printf("\nQ (mean queue lengths, indicative interval):\n");
		for (k = 0; k < M; k++) for (r = 0; r < R; r++) {
			char lbl[32]; snprintf(lbl, sizeof lbl, "Q[%d][%d]", k+1, r+1);
			print_bracket(lbl, om.Q[k*R+r], op.Q[k*R+r]);
		}
	}
	printf("=========================================\n");
	free(op.X); free(op.Q); free(om.X); free(om.Q);
	return 0;
}
