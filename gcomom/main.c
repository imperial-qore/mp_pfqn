/* gcomom - CoMoM basis with multi-branch queue removal (b = r-1, depth 1).
 *
 * While class c (0-based) is swept, classes s < c sit at their full
 * population N_s, class c at n, classes > c at 0.  A sub-model is a vector
 * v of station copies (v_k <= mi_k; mi is native, a station of multiplicity
 * m being m identical single-server queues).  For a decrement vector d over
 * the classes s < c with |d| <= D(v) it carries
 *     P[d]   = G_v(n - d)            Q[k,d] = G_{v+e_k}(n - d)
 * and advances them one population step of class c at a time with these
 * exact product-form relations:
 *   class-c PC  n P[d] = Z_c P'[d] + sum_k v_k L_kc Q'[k,d]          (explicit)
 *   CE          Q[k,d] - sum_{s<c} L_ks Q[k,d+1_s] = P[d] + L_kc Q'[k,d]
 *   PC  (s<c)   sum_k v_k L_ks Q[k,d+1_s] = (N_s-d_s) P[d] - Z_s P[d+1_s]
 *   removal of one copy of station j (k != j)
 *               Q[k,d] - sum_{s<c} L_js Q[k,d+1_s] = Q^{v-e_j}[k,d] + L_jc Q'[k,d]
 * (' = previous step).  CE/PC/removal rows are written for |d| <= D-1,
 * unknowns are Q[k,d] for |d| <= D.  Subtracting CE from removal gives
 * sum_s (L_ks - L_js) Q[k,d+1_s] on the right-hand side, so removing r-1
 * stations determines the depth-1 values of each station directly; the PC
 * rows supply the one row each removed station lacks.  Models whose support
 * has >= R-1 stations therefore run at depth D = 1, removing one copy of each
 * of their top R-1 support stations; the closure of that rule is the
 * family of sub-models.  Smaller supports cannot reach depth 1 and run at
 * CoMoM depth D = |support|, with the removal rows as extra rows; a child
 * then always has depth >= parent depth - 1, which is what its removal rows
 * read.  Classes idle at every support station add the exact rows of their
 * decoupled population constraint.
 *
 * The coefficient matrix of a sub-model depends on (L, v, Z, N) but not on
 * the step, so once per class a full set of independent rows is chosen
 * (modulo a 61-bit prime, exactly over Q if that is inconclusive), the
 * square system is split into its independent blocks and each block is
 * LU-factored exactly; a step is then one back-substitution per block.
 * A rank-deficient sub-model makes the solve report "singular", and
 * util/singular.h takes over (class orders, folded stations, aggregated
 * classes, exact interpolation of a perturbation).  Every result is exact. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include "gmpla.h"
#include "util.h"
#include "singular.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static bool STATS = false;
static bool NO_FALLBACK = false;   /* -X: no convolution sub-models, no singular.h tiers */

static int NTHREADS = 1;            /* -j: sub-models of one level in parallel */

/* wall clock: with -j the CPU time of the process sums the threads */
static double cputime(void)
{
	struct timeval tv; gettimeofday(&tv, NULL);
	return tv.tv_sec + 1e-6*tv.tv_usec;
}

static void* xcalloc(size_t n, size_t s) { void* p = calloc(n ? n : 1, s); if (!p) { fprintf(stderr, "gcomom: out of memory\n"); exit(2); } return p; }

/* ------------------------------------------------ decrement vectors */
/* dl[c] lists the vectors d over c classes with |d| <= Dmax, grouped by |d|,
 * so the vectors with |d| <= D are a prefix of length ndc(c,D) and an index
 * is the same at every depth.  up[c][i*c+s] = index of d_i + 1_s, or -1. */
typedef struct { int c, n; int* d; int* up; int* tot; } declist;

static void dl_rec(declist* L, int c, int pos, int left, int* cur)
{
	if (pos == c - 1 || c == 0) {
		if (c > 0) cur[c-1] = left;
		memcpy(L->d + (size_t)L->n*c, cur, c*sizeof(int));
		L->n++;
		return;
	}
	for (int x = left; x >= 0; x--) { cur[pos] = x; dl_rec(L, c, pos+1, left-x, cur); }
}

static long ndc(int c, int D) { return nck(D + c, c); }

static int dl_find(const declist* L, const int* d, int len)
{
	for (int i = 0; i < len; i++) if (!memcmp(L->d + (size_t)i*L->c, d, L->c*sizeof(int))) return i;
	return -1;
}

static declist* dl_build(int c, int Dmax)
{
	declist* L = xcalloc(1, sizeof(declist));
	long n = ndc(c, Dmax);
	L->c = c; L->d = xcalloc((size_t)n*(c ? c : 1), sizeof(int));
	L->tot = xcalloc(n, sizeof(int)); L->up = xcalloc((size_t)n*(c ? c : 1), sizeof(int));
	int* cur = xcalloc(c ? c : 1, sizeof(int));
	if (c == 0) L->n = 1;
	else for (int t = 0; t <= Dmax; t++) dl_rec(L, c, 0, t, cur);
	for (int i = 0; i < L->n; i++) {
		int t = 0; for (int s = 0; s < c; s++) t += L->d[(size_t)i*c+s];
		L->tot[i] = t;
		for (int s = 0; s < c; s++) {
			if (t + 1 > Dmax) { L->up[(size_t)i*c+s] = -1; continue; }
			memcpy(cur, L->d + (size_t)i*c, c*sizeof(int)); cur[s]++;
			L->up[(size_t)i*c+s] = dl_find(L, cur, n);
		}
	}
	free(cur);
	return L;
}

/* ------------------------------------------------ sub-models */
typedef struct { mpq_t* P; mpq_t* Q; } table;   /* P[nd], Q[nt*nd] */

enum { RCE, RPC, RHOM, RRM };
typedef struct { int kind, a, b, i; } rowdesc;  /* CE: a=ti; PC: a=s; RM: a=removal index, b=ti */

typedef struct {
	int n;           /* unknowns of the block */
	int* unk;        /* their indices */
	rowdesc* rows;   /* the selected rows, one per unknown */
	mpq_mat_t LU; int* idx;
	mpq_vec_t b;
} block;

typedef struct {
	int* v; int nt; int* T; int* pos; int D, tot;
	int nB; int* Bj; int* child;
	/* per class */
	int nd, nunk; int* uid;
	int nblk; block* blk;
	table** hist; int hlen;   /* hist[0] newest; at most D+1 */
	table* cur;
	/* singular sub-model with an affordable population lattice: its values
	 * come from exact convolution instead (lat = G_v scaled by prod n_s!,
	 * latk[ti] = G_{v+e_k} likewise) */
	int conv; mpz_t* lat; mpz_t* latk;
} gmodel;

typedef struct {
	int M, R; int* N; mpz_t* Z; mpz_t** L; int* mi;
	/* fractional multiplicity e_k = mi[k]/mden (util.h, qn_mult) splits into a
	 * base fb[k]/mden in [0,1), never removed, plus v_k removable copies, the
	 * integers of the copy vectors; fb = 0 and mden = 1 for integer models */
	int mden; int* fb; int* vfull;
	gmodel* m; int nm, full;
	declist** dl; int Dmax;
	long unk_total;
	long latn; long* rad;     /* population lattice size and mixed-radix strides */
	int nconv;
} gctx;

/* largest lattice-points x (support+1) a convolution sub-model may use */
#define CONV_MAX 4000000L

static table* tb_new(int nd, int nt)
{
	table* t = xcalloc(1, sizeof(table));
	t->P = xcalloc(nd, sizeof(mpq_t)); t->Q = xcalloc((size_t)nd*(nt ? nt : 1), sizeof(mpq_t));
	for (int i = 0; i < nd; i++) mpq_init(t->P[i]);
	for (long i = 0; i < (long)nd*nt; i++) mpq_init(t->Q[i]);
	return t;
}
static void tb_free(table* t, int nd, int nt)
{
	if (!t) return;
	for (int i = 0; i < nd; i++) mpq_clear(t->P[i]);
	for (long i = 0; i < (long)nd*nt; i++) mpq_clear(t->Q[i]);
	free(t->P); free(t->Q); free(t);
}

/* family: open-addressing hash on copy vectors */
typedef struct { int cap; int* key; } vhash;
static uint64_t vh_code(const int* v, int M) { uint64_t h = 1469598103934665603ULL; for (int k = 0; k < M; k++) { h ^= (uint64_t)(v[k]+1); h *= 1099511628211ULL; } return h; }

static void vh_insert(vhash* H, gctx* g, int id)
{
	uint64_t h = vh_code(g->m[id].v, g->M) & (H->cap-1);
	while (H->key[h] >= 0) h = (h+1) & (H->cap-1);
	H->key[h] = id;
}

/* index of copy vector v in the family, added if new */
static int vh_find_or_add(vhash* H, gctx* g, int* cap, const int* v)
{
	int M = g->M;
	uint64_t h = vh_code(v, M) & (H->cap-1);
	while (H->key[h] >= 0) {
		if (!memcmp(g->m[H->key[h]].v, v, M*sizeof(int))) return H->key[h];
		h = (h+1) & (H->cap-1);
	}
	if (g->nm == *cap) {
		*cap *= 2;
		g->m = realloc(g->m, (size_t)*cap*sizeof(gmodel));
		if (!g->m) { fprintf(stderr, "gcomom: out of memory\n"); exit(2); }
		memset(g->m + g->nm, 0, (size_t)(*cap - g->nm)*sizeof(gmodel));
	}
	int id = g->nm++;
	g->m[id].v = xcalloc(M, sizeof(int)); memcpy(g->m[id].v, v, M*sizeof(int));
	if (2*g->nm > H->cap) {                /* grow and rehash */
		free(H->key); H->cap *= 2; H->key = xcalloc(H->cap, sizeof(int));
		for (int z = 0; z < H->cap; z++) H->key[z] = -1;
		for (int z = 0; z < g->nm; z++) vh_insert(H, g, z);
	} else H->key[h] = id;
	return id;
}

/* the full model and the closure of: remove one copy of each of the top
 * min(R-1, |support|) support stations */
static void fam_build(gctx* g)
{
	int M = g->M, R = g->R, cap = 1024;
	g->m = xcalloc(cap, sizeof(gmodel)); g->nm = 0;
	vhash H = { 4096, xcalloc(4096, sizeof(int)) };
	for (int z = 0; z < H.cap; z++) H.key[z] = -1;
	int* tmp = xcalloc(M, sizeof(int));
	g->full = vh_find_or_add(&H, g, &cap, g->vfull);
	for (int id = 0; id < g->nm; id++) {
		gmodel* x = &g->m[id];
		int nt = 0;
		for (int k = 0; k < M; k++) if (x->v[k] > 0 || g->fb[k] > 0) nt++;
		x->nt = nt; x->T = xcalloc(nt, sizeof(int)); x->pos = xcalloc(M, sizeof(int));
		nt = 0; x->tot = 0;
		int nrem = 0;   /* support stations with a removable copy (all of them in integer models) */
		for (int k = 0; k < M; k++) { x->pos[k] = -1; if (x->v[k] > 0 || g->fb[k] > 0) { x->pos[k] = nt; x->T[nt++] = k; x->tot += x->v[k]; if (x->v[k] > 0) nrem++; } }
		x->D = (nrem >= R-1) ? 1 : (nt > 0 ? nt : 1);
		x->nB = MIN(R-1, nrem);
		x->Bj = xcalloc(x->nB, sizeof(int)); x->child = xcalloc(x->nB, sizeof(int));
		int* rem = xcalloc(nrem ? nrem : 1, sizeof(int));
		{ int z = 0; for (int ti = 0; ti < nt; ti++) if (x->v[x->T[ti]] > 0) rem[z++] = x->T[ti]; }
		for (int b = 0; b < x->nB; b++) {
			int j = rem[nrem - x->nB + b], any = 0;
			memcpy(tmp, x->v, M*sizeof(int)); tmp[j]--;
			for (int k = 0; k < M; k++) if (tmp[k] || g->fb[k]) { any = 1; break; }
			x->Bj[b] = j;
			/* vh_find_or_add may move g->m: re-fetch x after it */
			int ch = any ? vh_find_or_add(&H, g, &cap, tmp) : -1;   /* the empty model is never read */
			x = &g->m[id];
			x->child[b] = ch;
		}
		free(rem);
	}
	free(tmp); free(H.key);
	g->Dmax = 1; for (int i = 0; i < g->nm; i++) g->Dmax = MAX(g->Dmax, g->m[i].D);
}

/* order by total copies so that children precede parents */
static int* fam_order(gctx* g)
{
	int* o = xcalloc(g->nm, sizeof(int)), n = 0, maxt = 0;
	for (int i = 0; i < g->nm; i++) maxt = MAX(maxt, g->m[i].tot);
	for (int t = 0; t <= maxt; t++) for (int i = 0; i < g->nm; i++) if (g->m[i].tot == t) o[n++] = i;
	return o;
}

/* ------------------------------------------------ row selection */
#define P61 2305843009213693951ULL
#define P62 4611686018427387847ULL
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)((unsigned __int128)a*b % p); }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t p) { uint64_t r = 1; while (e) { if (e&1) r = mulmod(r,a,p); a = mulmod(a,a,p); e >>= 1; } return r; }

/* Choose rows of A (nr x nc, integers) spanning its row space; returns the
 * rank, sel[] the chosen rows.  Modulo p first: full rank mod p proves full
 * rank over Q for those rows.  Otherwise retry with a second prime and
 * finally decide exactly over Q. */
static int select_rows_modp(mpz_t* A, int nr, int nc, int* sel, uint64_t p)
{
	uint64_t* bas = xcalloc((size_t)nc*nc, sizeof(uint64_t));
	int* pc = xcalloc(nc, sizeof(int)), rk = 0;
	uint64_t* x = xcalloc(nc, sizeof(uint64_t));
	for (int r = 0; r < nr && rk < nc; r++) {
		for (int j = 0; j < nc; j++) x[j] = mpz_fdiv_ui(A[(size_t)r*nc+j], p);
		for (int b = 0; b < rk; b++) {
			uint64_t f = x[pc[b]]; if (!f) continue;
			uint64_t* br = bas + (size_t)b*nc;
			for (int j = 0; j < nc; j++) if (br[j]) x[j] = (x[j] + p - mulmod(f, br[j], p)) % p;
		}
		int c = -1; for (int j = 0; j < nc; j++) if (x[j]) { c = j; break; }
		if (c < 0) continue;
		uint64_t inv = powmod(x[c], p-2, p);
		uint64_t* br = bas + (size_t)rk*nc;
		for (int j = 0; j < nc; j++) br[j] = mulmod(x[j], inv, p);
		pc[rk] = c; sel[rk++] = r;
	}
	free(bas); free(pc); free(x);
	return rk;
}

static int select_rows_exact(mpz_t* A, int nr, int nc, int* sel)
{
	mpq_t* bas = xcalloc((size_t)nc*nc, sizeof(mpq_t)); mpq_t* x = xcalloc(nc, sizeof(mpq_t));
	for (long i = 0; i < (long)nc*nc; i++) mpq_init(bas[i]);
	for (int j = 0; j < nc; j++) mpq_init(x[j]);
	int* pc = xcalloc(nc, sizeof(int)), rk = 0;
	mpq_t f, t; mpq_init(f); mpq_init(t);
	for (int r = 0; r < nr && rk < nc; r++) {
		for (int j = 0; j < nc; j++) mpq_set_z(x[j], A[(size_t)r*nc+j]);
		for (int b = 0; b < rk; b++) {
			if (!mpq_sgn(x[pc[b]])) continue;
			mpq_set(f, x[pc[b]]);
			for (int j = 0; j < nc; j++) if (mpq_sgn(bas[(size_t)b*nc+j])) { mpq_mul(t, f, bas[(size_t)b*nc+j]); mpq_sub(x[j], x[j], t); }
		}
		int c = -1; for (int j = 0; j < nc; j++) if (mpq_sgn(x[j])) { c = j; break; }
		if (c < 0) continue;
		mpq_inv(f, x[c]);
		for (int j = 0; j < nc; j++) mpq_mul(bas[(size_t)rk*nc+j], x[j], f);
		pc[rk] = c; sel[rk++] = r;
	}
	for (long i = 0; i < (long)nc*nc; i++) mpq_clear(bas[i]);
	for (int j = 0; j < nc; j++) mpq_clear(x[j]);
	free(bas); free(x); free(pc); mpq_clear(f); mpq_clear(t);
	return rk;
}

static int select_rows(mpz_t* A, int nr, int nc, int* sel)
{
	if (select_rows_modp(A, nr, nc, sel, P61) == nc) return nc;
	if (select_rows_modp(A, nr, nc, sel, P62) == nc) return nc;
	return select_rows_exact(A, nr, nc, sel);
}

/* ------------------------------------------------ per-class setup */
static int uf_find(int* p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

static void blocks_free(gmodel* x)
{
	for (int b = 0; b < x->nblk; b++) {
		block* B = &x->blk[b];
		for (int i = 0; i < B->n; i++) { for (int j = 0; j < B->n; j++) mpq_clear(B->LU[i][j]); free(B->LU[i]); }
		free(B->LU); free(B->idx); free(B->unk); free(B->rows);
		for (int i = 0; i <= B->n; i++) mpq_clear(B->b[i]);
		free(B->b);
	}
	free(x->blk); x->blk = NULL; x->nblk = 0;
}

static bool valid_d(gctx* g, int c, int i)
{
	const int* d = g->dl[c]->d + (size_t)i*c;
	for (int s = 0; s < c; s++) if (d[s] > g->N[s]) return false;
	return true;
}

/* Build the rows of model x for class c and factor its blocks.
 * Returns 0, or 1 if the coefficient matrix is rank deficient. */
static int setup_class(gctx* g, gmodel* x, int c)
{
	declist* dl = g->dl[c];
	int nd = x->nd, nt = x->nt;   /* set by the caller: nck() caches into a global table, not thread safe */
	free(x->uid); x->uid = xcalloc((size_t)nt*nd, sizeof(int));
	int nunk = 0;
	for (int ti = 0; ti < nt; ti++) for (int i = 0; i < nd; i++) x->uid[ti*nd+i] = valid_d(g, c, i) ? nunk++ : -1;
	x->nunk = nunk;
	#pragma omp atomic
	g->unk_total += nunk;
	/* rows */
	int maxrows = 0, nr = 0;
	for (int i = 0; i < nd; i++) if (dl->tot[i] <= x->D-1) maxrows += nt + c + 2*c*nt + x->nB*nt;
	rowdesc* rd = xcalloc(maxrows, sizeof(rowdesc));
	mpz_t* A = xcalloc((size_t)maxrows*(nunk ? nunk : 1), sizeof(mpz_t));
	for (long z = 0; z < (long)maxrows*nunk; z++) mpz_init(A[z]);
	mpz_t w; mpz_init(w);
#define UNK(ti, ii) (((ii) >= 0 && (ii) < nd) ? x->uid[(ti)*nd+(ii)] : -1)
#define ADD(col, val) do { int cc_ = (col); if (cc_ >= 0) mpz_add(A[(size_t)nr*nunk+cc_], A[(size_t)nr*nunk+cc_], (val)); } while (0)
	for (int i = 0; i < nd; i++) {
		if (dl->tot[i] > x->D-1 || !valid_d(g, c, i)) continue;
		const int* d = dl->d + (size_t)i*c;
		for (int ti = 0; ti < nt; ti++) {                     /* CE */
			int k = x->T[ti];
			mpz_set_ui(w, 1); ADD(UNK(ti, i), w);
			for (int s = 0; s < c; s++) { mpz_neg(w, g->L[k][s]); ADD(UNK(ti, dl->up[(size_t)i*c+s]), w); }
			rd[nr++] = (rowdesc){ RCE, ti, 0, i };
		}
		for (int s = 0; s < c; s++) {                          /* PC */
			int u = dl->up[(size_t)i*c+s];
			/* (fb_k/mden + v_k) L_ks, the row scaled by mden (its right-hand side too, in step_model) */
			for (int ti = 0; ti < nt; ti++) { int k = x->T[ti]; mpz_mul_si(w, g->L[k][s], g->fb[k] + (long)g->mden * x->v[k]); ADD(UNK(ti, u), w); }
			rd[nr++] = (rowdesc){ RPC, s, 0, i };
		}
		for (int s = 0; s < c; s++) {                          /* class s idle in the support */
			bool idle = true;
			for (int ti = 0; ti < nt; ti++) if (mpz_sgn(g->L[x->T[ti]][s])) { idle = false; break; }
			if (!idle) continue;
			int u = dl->up[(size_t)i*c+s];
			for (int ti = 0; ti < nt; ti++) {
				if (mpz_sgn(g->Z[s]) == 0) {
					/* no demand and no think time: G = 0 while class s has jobs */
					int dd[2] = { i, u };
					for (int h = 0; h < 2; h++) {
						int ii = dd[h];
						if (ii < 0 || ii >= nd || !valid_d(g, c, ii)) continue;
						if (dl->d[(size_t)ii*c+s] >= g->N[s]) continue;
						mpz_set_ui(w, 1); ADD(UNK(ti, ii), w);
						rd[nr++] = (rowdesc){ RHOM, 0, 0, i };
					}
				} else {
					/* Z_s Q[k,d+1_s] = (N_s - d_s) Q[k,d] */
					ADD(UNK(ti, u), g->Z[s]);
					mpz_set_si(w, -(g->N[s] - d[s])); ADD(UNK(ti, i), w);
					rd[nr++] = (rowdesc){ RHOM, 0, 0, i };
				}
			}
		}
		for (int b = 0; b < x->nB; b++) {                      /* queue removal */
			int j = x->Bj[b];
			if (x->child[b] < 0) continue;
			for (int ti = 0; ti < nt; ti++) {
				if (x->T[ti] == j) continue;
				mpz_set_ui(w, 1); ADD(UNK(ti, i), w);
				for (int s = 0; s < c; s++) { mpz_neg(w, g->L[j][s]); ADD(UNK(ti, dl->up[(size_t)i*c+s]), w); }
				rd[nr++] = (rowdesc){ RRM, b, ti, i };
			}
		}
	}
#undef ADD
#undef UNK
	int* sel = xcalloc(nunk ? nunk : 1, sizeof(int));
	int rk = select_rows(A, nr, nunk, sel);
	int rc = 0;
	if (rk < nunk) { rc = 1; goto done; }
	/* blocks: unknowns linked by a selected row */
	{
		int* par = xcalloc(nunk, sizeof(int));
		for (int u = 0; u < nunk; u++) par[u] = u;
		for (int q = 0; q < nunk; q++) {
			int first = -1;
			for (int u = 0; u < nunk; u++) if (mpz_sgn(A[(size_t)sel[q]*nunk+u])) {
				if (first < 0) first = u; else par[uf_find(par, u)] = uf_find(par, first);
			}
		}
		int* bid = xcalloc(nunk, sizeof(int)), nb = 0;
		for (int u = 0; u < nunk; u++) bid[u] = -1;
		int* root2b = xcalloc(nunk, sizeof(int));
		for (int u = 0; u < nunk; u++) root2b[u] = -1;
		for (int u = 0; u < nunk; u++) { int r = uf_find(par, u); if (root2b[r] < 0) root2b[r] = nb++; bid[u] = root2b[r]; }
		blocks_free(x);
		x->nblk = nb; x->blk = xcalloc(nb, sizeof(block));
		int* cnt = xcalloc(nb, sizeof(int));
		for (int u = 0; u < nunk; u++) cnt[bid[u]]++;
		for (int b = 0; b < nb; b++) {
			x->blk[b].unk = xcalloc(cnt[b], sizeof(int)); x->blk[b].rows = xcalloc(cnt[b], sizeof(rowdesc));
		}
		int* lpos = xcalloc(nunk, sizeof(int));
		for (int u = 0; u < nunk; u++) { block* B = &x->blk[bid[u]]; lpos[u] = B->n; B->unk[B->n++] = u; }
		int* rcount = xcalloc(nb, sizeof(int));
		for (int b = 0; b < nb; b++) {
			block* B = &x->blk[b];
			B->LU = mpq_matzeros(B->n, B->n);
		}
		for (int q = 0; q < nunk; q++) {
			int r = sel[q], b = -1;
			for (int u = 0; u < nunk; u++) if (mpz_sgn(A[(size_t)r*nunk+u])) { b = bid[u]; break; }
			block* B = &x->blk[b];
			int rr = rcount[b]++;
			B->rows[rr] = rd[r];
			for (int u = 0; u < nunk; u++) if (mpz_sgn(A[(size_t)r*nunk+u])) mpq_set_z(B->LU[rr][lpos[u]], A[(size_t)r*nunk+u]);
		}
		for (int b = 0; b < nb && !rc; b++) {
			block* B = &x->blk[b];
			if (rcount[b] != B->n) { rc = 1; break; }   /* cannot happen for a nonsingular selection */
			B->idx = mpq_ludcmp(B->LU, B->n);
			if (!B->idx) rc = 1;
			B->b = mpq_vec(B->n, 0, 1);
		}
		free(par); free(bid); free(root2b); free(cnt); free(lpos); free(rcount);
	}
done:
	for (long z = 0; z < (long)maxrows*nunk; z++) mpz_clear(A[z]);
	free(A); free(rd); free(sel); mpz_clear(w);
	return rc;
}

/* ------------------------------------------------ convolution sub-models */
/* Exact Buzen convolution of sub-model x over the whole population lattice,
 * in integers: with h(n) = G(n) prod_s n_s!, the delay gives prod_s Z_s^{n_s}
 * and each single-server queue copy of station k
 *     h(n) = h_prev(n) + sum_s L_ks n_s h(n - 1_s).
 * lat holds G_v, latk[ti] holds G_{v+e_k} for the support stations. */
static void conv_one(gctx* g, mpz_t* h, int k)
{
	int R = g->R; mpz_t t; mpz_init(t);
	int* p = xcalloc(R, sizeof(int));
	for (long i = 0; i < g->latn; i++) {   /* lexicographic order: n - 1_s precedes n */
		long rem = i;
		for (int s = R-1; s >= 0; s--) { p[s] = rem / g->rad[s]; rem %= g->rad[s]; }
		for (int s = 0; s < R; s++) if (p[s] > 0) {
			mpz_mul_ui(t, g->L[k][s], p[s]); mpz_addmul(h[i], t, h[i - g->rad[s]]);
		}
	}
	free(p); mpz_clear(t);
}

static void conv_build(gctx* g, gmodel* x)
{
	int R = g->R; long n = g->latn;
	x->lat = xcalloc(n, sizeof(mpz_t)); x->latk = xcalloc((size_t)n*(x->nt ? x->nt : 1), sizeof(mpz_t));
	int* p = xcalloc(R, sizeof(int)); mpz_t t; mpz_init(t);
	for (long i = 0; i < n; i++) {                        /* delay */
		long rem = i; mpz_init_set_ui(x->lat[i], 1);
		for (int s = R-1; s >= 0; s--) { p[s] = rem / g->rad[s]; rem %= g->rad[s]; }
		for (int s = 0; s < R; s++) { mpz_pow_ui(t, g->Z[s], p[s]); mpz_mul(x->lat[i], x->lat[i], t); }
	}
	for (int ti = 0; ti < x->nt; ti++)
		for (int cpy = 0; cpy < x->v[x->T[ti]]; cpy++) conv_one(g, x->lat, x->T[ti]);
	for (int ti = 0; ti < x->nt; ti++) {
		mpz_t* hk = x->latk + (size_t)ti*n;
		for (long i = 0; i < n; i++) mpz_init_set(hk[i], x->lat[i]);
		conv_one(g, hk, x->T[ti]);
	}
	free(p); mpz_clear(t);
	x->conv = 1;   /* the caller counts g->nconv */
}

static void conv_free(gctx* g, gmodel* x)
{
	if (!x->conv) return;
	for (long i = 0; i < g->latn; i++) mpz_clear(x->lat[i]);
	for (long i = 0; i < (long)g->latn*x->nt; i++) mpz_clear(x->latk[i]);
	free(x->lat); free(x->latk);
}

/* table of a convolution sub-model at class c, class-c population n */
static table* conv_table(gctx* g, gmodel* x, int c, int n)
{
	declist* dl = g->dl[c];
	int nd = x->nd, nt = x->nt;
	table* T = tb_new(nd, nt);
	mpz_t f; mpz_init(f);
	for (int i = 0; i < nd; i++) {
		if (!valid_d(g, c, i)) continue;
		long idx = (long)n * g->rad[c];
		mpz_set_ui(f, 1);
		for (int s = 0; s < c; s++) {
			int ps = g->N[s] - dl->d[(size_t)i*c+s];
			idx += (long)ps * g->rad[s];
			mpz_t fs; mpz_init(fs); mpz_fac_ui(fs, ps); mpz_mul(f, f, fs); mpz_clear(fs);
		}
		{ mpz_t fs; mpz_init(fs); mpz_fac_ui(fs, n); mpz_mul(f, f, fs); mpz_clear(fs); }
		mpq_set_num(T->P[i], x->lat[idx]); mpq_set_den(T->P[i], f); mpq_canonicalize(T->P[i]);
		for (int ti = 0; ti < nt; ti++) {
			mpq_ptr q = T->Q[ti*nd+i];
			mpq_set_num(q, x->latk[(size_t)ti*g->latn + idx]); mpq_set_den(q, f); mpq_canonicalize(q);
		}
	}
	mpz_clear(f);
	return T;
}

/* ------------------------------------------------ one step */
static void step_model(gctx* g, gmodel* x, int c, int n)
{
	if (x->conv) { x->cur = conv_table(g, x, c, n); return; }
	declist* dl = g->dl[c];
	int nd = x->nd, nt = x->nt;
	table* pv = x->hist[0];
	table* T = tb_new(nd, nt);
	mpq_t t, w; mpq_init(t); mpq_init(w);
	/* plain constants from the class-c population constraint */
	for (int i = 0; i < nd; i++) {
		if (!valid_d(g, c, i)) continue;
		mpq_set_z(w, g->Z[c]); mpq_mul(T->P[i], w, pv->P[i]);
		for (int ti = 0; ti < nt; ti++) {
			int k = x->T[ti];
			mpz_mul_si(mpq_numref(w), g->L[k][c], g->fb[k] + (long)g->mden * x->v[k]); mpz_set_ui(mpq_denref(w), g->mden); mpq_canonicalize(w);
			mpq_mul(t, w, pv->Q[ti*nd+i]); mpq_add(T->P[i], T->P[i], t);
		}
		mpq_set_si(w, 1, n); mpq_mul(T->P[i], T->P[i], w);
	}
	/* blocks */
	for (int b = 0; b < x->nblk; b++) {
		block* B = &x->blk[b];
		for (int r = 0; r < B->n; r++) {
			rowdesc* R = &B->rows[r];
			mpq_ptr y = B->b[r];
			int i = R->i;
			switch (R->kind) {
			case RCE: {
				int k = x->T[R->a];
				mpq_set_z(w, g->L[k][c]); mpq_mul(y, w, pv->Q[R->a*nd+i]); mpq_add(y, y, T->P[i]);
				break; }
			case RPC: {
				int s = R->a, u = dl->up[(size_t)i*c+s];
				mpq_set_si(w, g->N[s] - dl->d[(size_t)i*c+s], 1); mpq_mul(y, w, T->P[i]);
				if (u >= 0 && u < nd) { mpq_set_z(w, g->Z[s]); mpq_mul(t, w, T->P[u]); mpq_sub(y, y, t); }
				if (g->mden > 1) { mpq_set_si(w, g->mden, 1); mpq_mul(y, y, w); }   /* row scaled by mden */
				break; }
			case RHOM:
				mpq_set_ui(y, 0, 1);
				break;
			case RRM: {
				int j = x->Bj[R->a], k = x->T[R->b];
				gmodel* ch = &g->m[x->child[R->a]];
				mpq_set_z(w, g->L[j][c]); mpq_mul(y, w, pv->Q[R->b*nd+i]);
				mpq_add(y, y, ch->cur->Q[ch->pos[k]*ch->nd + i]);
				break; }
			}
		}
		mpq_lubksb(B->LU, B->b, B->n, B->idx);
	}
	/* scatter the block solutions back to Q */
	{
		int* slotb = xcalloc(x->nunk ? x->nunk : 1, sizeof(int)), *slotr = xcalloc(x->nunk ? x->nunk : 1, sizeof(int));
		for (int b = 0; b < x->nblk; b++) for (int r = 0; r < x->blk[b].n; r++) { slotb[x->blk[b].unk[r]] = b; slotr[x->blk[b].unk[r]] = r; }
		for (int ti = 0; ti < nt; ti++) for (int i = 0; i < nd; i++) {
			int u = x->uid[ti*nd+i];
			if (u >= 0) mpq_set(T->Q[ti*nd+i], x->blk[slotb[u]].b[slotr[u]]);
		}
		free(slotb); free(slotr);
	}
	mpq_clear(t); mpq_clear(w);
	x->cur = T;
}

static void push_hist(gmodel* x, table* T)
{
	int cap = x->D + 1;
	if (x->hlen == cap) { tb_free(x->hist[cap-1], x->nd, x->nt); x->hlen--; }
	memmove(x->hist + 1, x->hist, x->hlen*sizeof(table*));
	x->hist[0] = T; x->hlen++;
}

static void hist_clear(gmodel* x, int nd)
{
	for (int h = 0; h < x->hlen; h++) tb_free(x->hist[h], nd, x->nt);
	x->hlen = 0;
}

/* ------------------------------------------------ the solver */
/* Exported constants of q (util/singular.h layout).  0 solved, 1 singular. */
static int gcomom_solve(qnmodel* q, mpq_t* out, void* ctx)
{
	(void)ctx;
	double t0 = cputime();
	gctx G; memset(&G, 0, sizeof G);
	gctx* g = &G;
	g->M = q->M; g->R = q->R; g->N = q->N; g->Z = q->Z; g->L = q->L; g->mi = q->mi;
	g->mden = qn_mden(q);
	g->fb = xcalloc(q->M, sizeof(int)); g->vfull = xcalloc(q->M, sizeof(int));
	for (int k = 0; k < q->M; k++) { g->fb[k] = q->mi[k] % g->mden; g->vfull[k] = q->mi[k] / g->mden; }
	int M = g->M, R = g->R, rc = 0;
	int Ntot = 0; for (int r = 0; r < R; r++) Ntot += g->N[r];
	nckinit(MAX(M + R + 8, 2*R + 8), R + 2);
	bool anyst = false; for (int k = 0; k < M; k++) if (g->mi[k] > 0) anyst = true;
	if (!anyst) { fprintf(stderr, "gcomom: the model has no station\n"); return 2; }
	fam_build(g);
	int* ord = fam_order(g);
	g->dl = xcalloc(R, sizeof(declist*));
	for (int c = 0; c < R; c++) g->dl[c] = dl_build(c, g->Dmax);
	/* population lattice, for sub-models solved by convolution; latn = 0 if huge */
	g->rad = xcalloc(R, sizeof(long)); g->latn = 1;
	for (int r = 0; r < R; r++) {
		g->rad[r] = g->latn;
		if (g->latn > CONV_MAX / (g->N[r] + 1)) { g->latn = 0; break; }
		g->latn *= g->N[r] + 1;
	}
	/* level boundaries in ord (fam_order sorts by total copies) */
	int* lev = xcalloc(g->nm + 1, sizeof(int)), nlev = 0;
	for (int oi = 0; oi < g->nm; oi++) if (oi == 0 || g->m[ord[oi]].tot != g->m[ord[oi-1]].tot) lev[nlev++] = oi;
	lev[nlev] = g->nm;
	int* sing = xcalloc(g->nm, sizeof(int)), nsing = 0;
	/* class 0, n = 0: zero population */
	for (int i = 0; i < g->nm; i++) {
		gmodel* x = &g->m[i];
		x->hist = xcalloc(x->D + 1, sizeof(table*));
		x->nd = 1;
		table* T = tb_new(1, x->nt);
		mpq_set_ui(T->P[0], 1, 1);
		for (int ti = 0; ti < x->nt; ti++) mpq_set_ui(T->Q[ti], 1, 1);
		push_hist(x, T);
	}
	for (int c = 0; c < R && !rc; c++) {
		if (c > 0) {
			/* tables at n_c = 0 from the last steps of class c-1 */
			declist* dn = g->dl[c]; int Nl = g->N[c-1];
			for (int i = 0; i < g->nm; i++) {
				gmodel* x = &g->m[i];
				int ndo = x->nd, ndn = (int)ndc(c, x->D), nt = x->nt;
				table* T = tb_new(ndn, nt);
				for (int ii = 0; ii < ndn; ii++) {
					const int* d = dn->d + (size_t)ii*c;
					int tt = d[c-1];
					if (tt > Nl || !valid_d(g, c, ii)) continue;
					int src = dl_find(g->dl[c-1], d, ndo);
					if (src < 0 || tt >= x->hlen) continue;   /* |d'| <= D - t is always stored */
					table* S = x->hist[tt];
					mpq_set(T->P[ii], S->P[src]);
					for (int ti = 0; ti < nt; ti++) mpq_set(T->Q[ti*ndn+ii], S->Q[ti*ndo+src]);
				}
				hist_clear(x, ndo);
				x->nd = ndn;
				push_hist(x, T);
			}
		}
		if (g->N[c] == 0) continue;
		/* setup reads no other sub-model: all of them in parallel */
		#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
		for (int oi = 0; oi < g->nm; oi++) {
			gmodel* x = &g->m[ord[oi]];
			sing[oi] = x->conv ? 0 : setup_class(g, x, c);
		}
		nsing = 0;
		for (int oi = 0; oi < g->nm && !rc; oi++) {
			if (!sing[oi]) continue;
			gmodel* x = &g->m[ord[oi]];
			bool fits = !NO_FALLBACK && g->latn > 0 && g->latn <= CONV_MAX / (x->nt + 1);
			for (int ti = 0; ti < x->nt; ti++) if (g->fb[x->T[ti]]) fits = false;   /* conv_one adds whole copies only */
			if (STATS) {
				fprintf(stderr, "gcomom: singular sub-model v=(");
				for (int k = 0; k < M; k++) fprintf(stderr, "%d%s", x->v[k], k < M-1 ? "," : "");
				fprintf(stderr, ") at class %d%s\n", c+1, fits ? ", solved by convolution from here on" : "");
			}
			if (fits) { blocks_free(x); sing[nsing++] = ord[oi]; } else rc = 1;
		}
		if (rc) break;
		#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
		for (int z = 0; z < nsing; z++) conv_build(g, &g->m[sing[z]]);
		g->nconv += nsing;
		/* a step reads the children's new tables, which have fewer copies:
		 * the sub-models of one level (equal total copies) in parallel */
		for (int n = 1; n <= g->N[c] && !rc; n++) {
			for (int lv = 0; lv < nlev; lv++) {
				#pragma omp parallel for schedule(dynamic,1) num_threads(NTHREADS)
				for (int oi = lev[lv]; oi < lev[lv+1]; oi++) step_model(g, &g->m[ord[oi]], c, n);
			}
			for (int oi = 0; oi < g->nm; oi++) { gmodel* x = &g->m[ord[oi]]; push_hist(x, x->cur); x->cur = NULL; }
		}
	}
	if (!rc) {
		gmodel* f = &g->m[g->full];
		int c = R - 1, nd = f->nd;
		declist* dl = g->dl[c];
		table* T = f->hist[0];
		table* Pv = (g->N[c] > 0 && f->hlen > 1) ? f->hist[1] : NULL;
		mpq_set(out[0], T->P[0]);
		for (int s = 0; s < R; s++) {
			for (int k = 0; k < M; k++) mpq_set_ui(out[SG_IDX_GK(R,k,s)], 0, 1);
			mpq_set_ui(out[1+s], 0, 1);
			if (s < c) {
				int u = dl->up[s];
				mpq_set(out[1+s], T->P[u]);
				for (int k = 0; k < M; k++) if (f->pos[k] >= 0) mpq_set(out[SG_IDX_GK(R,k,s)], T->Q[f->pos[k]*nd+u]);
			} else if (Pv) {
				mpq_set(out[1+s], Pv->P[0]);
				for (int k = 0; k < M; k++) if (f->pos[k] >= 0) mpq_set(out[SG_IDX_GK(R,k,s)], Pv->Q[f->pos[k]*nd]);
			}
		}
		if (mpq_sgn(out[0]) == 0) rc = 1;
	}
	if (STATS)
		fprintf(stderr, "gcomom: %d sub-models (%d by convolution), depth <= %d, %ld unknowns set up, %.3f s%s\n",
		        g->nm, g->nconv, g->Dmax, g->unk_total, cputime() - t0, rc ? " (singular)" : "");
	/* free */
	for (int i = 0; i < g->nm; i++) {
		gmodel* x = &g->m[i];
		hist_clear(x, x->nd);
		conv_free(g, x);
		blocks_free(x);
		free(x->hist); free(x->v); free(x->T); free(x->pos); free(x->Bj); free(x->child); free(x->uid);
	}
	for (int c = 0; c < R; c++) { free(g->dl[c]->d); free(g->dl[c]->up); free(g->dl[c]->tot); free(g->dl[c]); }
	free(g->dl); free(g->m); free(ord); free(g->rad); free(lev); free(sing);
	free(g->fb); free(g->vfull);
	(void)Ntot;
	return rc;
}

/* ------------------------------------------------ output */
static void put_double(const mpq_t q, mpf_t f) { mpf_set_q(f, q); printf("%.15e", mpf_get_d(f)); }

/* X_r = G(N-1_r)/G and Q_kr = m_k L_kr G^{+k}(N-1_r)/G */
static void measure(qnmodel* q, mpq_t* ex, int k, int r, mpq_t out)
{
	int R = q->R;
	if (k < 0) { mpq_div(out, ex[1+r], ex[0]); return; }
	if (mpz_sgn(q->L[k][r]) == 0 || q->mi[k] == 0) { mpq_set_ui(out, 0, 1); return; }
	mpq_t w; mpq_init(w);
	mpq_div(out, ex[SG_IDX_GK(R,k,r)], ex[0]);
	mpq_set_z(w, q->L[k][r]); mpq_mul(out, out, w);
	qn_mult(w, q, k, 0); mpq_mul(out, out, w);
	mpq_clear(w);
}

static void usage(const char* p)
{
	printf("Usage: %s [options] model.qn\n", p);
	printf("  -l           log of the normalizing constant\n");
	printf("  -e           exact normalizing constant (numerator, denominator)\n");
	printf("  -g           normalizing constant as a double\n");
	printf("  -t           throughputs X_r\n");
	printf("  -q           mean queue lengths Q_kr, one row per station\n");
	printf("  -d           all metrics as exact rationals\n");
	printf("  -X, --no-perturb  fail on a singular system instead of the exact fallback\n");
	printf("  -j n         threads: sub-models of one level in parallel (default 1)\n");
	printf("  --stats      sub-model count, depth and timing on stderr\n");
	printf("  -h           this help\n");
}

int main(int argc, char** argv)
{
	bool lo = false, eo = false, go = false, to = false, qo = false, dbg = false, nop = false;
	char* file = NULL;
	for (int i = 1; i < argc; i++) {
		const char* a = argv[i];
		if (!strcmp(a, "-l")) lo = true;
		else if (!strcmp(a, "-e")) eo = true;
		else if (!strcmp(a, "-g")) go = true;
		else if (!strcmp(a, "-t")) to = true;
		else if (!strcmp(a, "-q")) qo = true;
		else if (!strcmp(a, "-d")) dbg = true;
		else if (!strcmp(a, "-X") || !strcmp(a, "--no-perturb")) { nop = true; NO_FALLBACK = true; }
		else if (!strcmp(a, "--stats")) STATS = true;
		else if (!strcmp(a, "-j") && i+1 < argc) { NTHREADS = atoi(argv[++i]); if (NTHREADS < 1) NTHREADS = 1; }
		else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0]); return 0; }
		else if (a[0] == '-') { fprintf(stderr, "gcomom: unknown option %s\n", a); usage(argv[0]); return 1; }
		else file = argv[i];
	}
	if (!file) { usage(argv[0]); return 1; }
	double t0 = cputime();
	qn_allow_rational_mi = 1;   /* fractional multiplicities (util.h, qn_mult) */
	qnmodel* q = readmodel(file);
	if (!q) { fprintf(stderr, "gcomom: cannot read %s\n", file); return 1; }
	qn_require_closed(q, "gcomom", 0);
	int R = q->R, M = q->M, nres = SG_NRES(q);
	/* a class with jobs, no think time and no demand anywhere has nowhere to
	 * go: G = 0 and every measure is undefined (the fallback would only
	 * interpolate its way to G = 0) */
	for (int r = 0; r < R; r++) {
		if (q->N[r] == 0 || mpz_sgn(q->Z[r])) continue;
		bool dem = false;
		for (int k = 0; k < M; k++) if (q->mi[k] > 0 && mpz_sgn(q->L[k][r])) { dem = true; break; }
		if (!dem) {
			fprintf(stderr, "gcomom: class %d has %d jobs but no think time and no demand at any station:\n"
			                "  G = 0 and the performance measures are undefined\n", r+1, q->N[r]);
			return 1;
		}
	}
	sg_qn_fix_nt(q);
	mpq_t* ex = sg_mpq_arr(nres);
	char tier[256] = "";
	int rc = gcomom_solve(q, ex, NULL);
	if (rc == 2) return 1;
	if (rc == 1) {
		if (nop) {
			fprintf(stderr, "gcomom: singular system and -X given; no exact answer\n");
			return 1;
		}
		sg_ops ops = { gcomom_solve, NULL, NULL, "" };
		rc = sg_solve_tiered(q, ex, tier, sizeof tier, &ops, 1);
		if (rc != 0) {
			fprintf(stderr, "gcomom: singular system, and every exact fallback failed (other class orders,\n"
			                "merged stations or classes, interpolation).  bin/ca solves any well-posed model.\n");
			return 1;
		}
		fprintf(stderr, "gcomom: singular system, solved exactly via %s\n", tier);
	}
	mpf_t f; mpf_init(f);
	mpq_t v, tot; mpq_init(v); mpq_init(tot);
	bool report = !lo && !eo && !go && !to && !qo;
	if (lo) printf("%.15e\n", mpq_log(ex[0]));
	else if (eo) gmp_printf("%Zd\n%Zd\n", mpq_numref(ex[0]), mpq_denref(ex[0]));
	else if (go) { put_double(ex[0], f); printf("\n"); }
	else if (to) for (int r = 0; r < R; r++) { measure(q, ex, -1, r, v); if (dbg) gmp_printf("%Qd\n", v); else { put_double(v, f); printf("\n"); } }
	else if (qo) for (int k = 0; k < M; k++) {
		for (int r = 0; r < R; r++) { measure(q, ex, k, r, v); if (dbg) gmp_printf("%Qd", v); else put_double(v, f); if (r < R-1) printf(" "); }
		printf("\n");
	}
	else if (report) {
		if (tier[0]) printf("\nSingular system, solved exactly via %s\n", tier);
		printf("\n========== Performance Metrics ==========\n");
		if (dbg) gmp_printf("G = %Qd\n", ex[0]);
		else { printf("G = "); put_double(ex[0], f); printf("\n"); printf("log(G) = %.15e\n", mpq_log(ex[0])); }
		printf("\nX (throughputs):\n");
		for (int r = 0; r < R; r++) {
			measure(q, ex, -1, r, v);
			if (dbg) gmp_printf("X[%d] = %Qd\n", r+1, v); else { printf("X[%d] = ", r+1); put_double(v, f); printf("\n"); }
		}
		printf("\nQ (mean queue lengths):\n");
		for (int k = 0; k < M; k++) {
			printf("Q[%d] =", k+1); mpq_set_ui(tot, 0, 1);
			for (int r = 0; r < R; r++) {
				measure(q, ex, k, r, v); mpq_add(tot, tot, v);
				if (dbg) gmp_printf("\t%Qd", v); else { printf("\t"); put_double(v, f); }
			}
			if (dbg) gmp_printf("\t(total: %Qd)\n", tot); else { printf("\t(total: "); put_double(tot, f); printf(")\n"); }
		}
		printf("=========================================\n");
		printf("Elapsed time (gCoMoM): %.6f s\n", cputime() - t0);
	}
	mpq_clear(v); mpq_clear(tot); mpf_clear(f);
	sg_mpq_arr_free(ex, nres);
	return 0;
}
