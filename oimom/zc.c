/* zc.c - zero-aware CoMoM for the derived networks of oimom.
 *
 * The class-oriented recursion of CoMoM (TSE'09): classes are added one at a
 * time; while class r goes from n_r = t-1 to t, the basis holds, for every
 * decrement c of the earlier classes E with |c| <= M,
 *   G_t(c)    = G(N_E - c, t)          X_t(k,c) = G^{+k}(N_E - c, t).
 * G_t comes from the class-r population constraint (PC) at the previous step,
 *   t G_t(c) = Z_r G_{t-1}(c) + sum_k m_k L_kr X_{t-1}(k,c),
 * and X_t solves the convolution equations (CE) and the PCs of the earlier
 * classes,
 *   X_t(k,c) - G_t(c) - sum_{s in E} L_ks X_t(k,c+e_s) = L_kr X_{t-1}(k,c),
 *   (N_s - c_s) G_t(c) - Z_s G_t(c+e_s) - sum_k m_k L_ks X_t(k,c+e_s) = 0,
 * written for |c| < M as in comom (the CE of a station with zero demand on
 * every earlier class is also written at |c| = M, where it closes).
 *
 * comom's square system is singular as soon as a station has a zero demand
 * on an earlier class: X_t(k,d) with |d| = M and d supported on the zero
 * classes of k appears in no equation (an "orphan").  Here
 *  - every unknown that no equation touches is dropped and marked unknown;
 *  - a G_t(c) that the PC above cannot give (it needs an orphan of the
 *    previous step) becomes an unknown of the system instead, determined by
 *    the PCs of the earlier classes through its Z_s term (and dropped too if
 *    every such Z_s is zero, since then nothing needs it);
 *  - independent rows are chosen by exact elimination (the system is
 *    consistent, possibly overdetermined), the chosen square block inverted
 *    once per class and reused at every step while the pattern of unknowns
 *    stays the same.
 * A value is only ever read where it is known, and the export reads a
 * dropped G^{+k}(N-e_r) only when L_kr = 0 (exported as 0: every caller
 * multiplies it by L_kr).  If a column finds no pivot, or a needed value is
 * unknown, the class order is changed; after every order fails the caller
 * falls back to CoMoM's tiered route.  Classes visiting no station (delay
 * only) are factored out first. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>

typedef struct {
	int nE, card, M;
	int* cv;     /* card x nE */
	int* nb;     /* card x nE: index of c+e_i, -1 outside |c|<=M */
	int* imp;    /* c_i > N_{E_i} for some i: every value is 0 */
	int* sum;    /* |c| */
} cset;

static int cs_find(const cset* C, const int* c)
{
	int i, j;
	for (i = 0; i < C->card; i++) {
		for (j = 0; j < C->nE; j++) if (C->cv[i * C->nE + j] != c[j]) break;
		if (j == C->nE) return i;
	}
	return -1;
}

static void cs_rec(cset* C, int* cur, int pos, int left, int* cap)
{
	int v;
	if (pos == C->nE) {
		if (C->card == *cap) { *cap *= 2; C->cv = (int*) realloc(C->cv, (size_t) *cap * (C->nE ? C->nE : 1) * sizeof(int)); }
		memcpy(C->cv + C->card * C->nE, cur, C->nE * sizeof(int));
		C->card++;
		return;
	}
	for (v = 0; v <= left; v++) { cur[pos] = v; cs_rec(C, cur, pos + 1, left - v, cap); }
}

static cset* cs_new(int nE, int M, const int* NE)
{
	cset* C = (cset*) calloc(1, sizeof(cset));
	int cap = 16, i, j, cur[32];
	C->nE = nE; C->M = M;
	C->cv = (int*) malloc((size_t) cap * (nE ? nE : 1) * sizeof(int));
	cs_rec(C, cur, 0, M, &cap);
	C->nb = (int*) malloc((size_t) C->card * (nE ? nE : 1) * sizeof(int));
	C->imp = (int*) calloc(C->card, sizeof(int));
	C->sum = (int*) calloc(C->card, sizeof(int));
	for (i = 0; i < C->card; i++) {
		int* c = C->cv + i * nE;
		for (j = 0; j < nE; j++) { C->sum[i] += c[j]; if (c[j] > NE[j]) C->imp[i] = 1; }
		for (j = 0; j < nE; j++) {
			if (C->sum[i] >= M) { C->nb[i * nE + j] = -1; continue; }
			memcpy(cur, c, nE * sizeof(int)); cur[j]++;
			C->nb[i * nE + j] = cs_find(C, cur);
		}
	}
	return C;
}

static void cs_free(cset* C) { free(C->cv); free(C->nb); free(C->imp); free(C->sum); free(C); }

/* values of one step: G[card], X[M*card] (X[k*card+c]), known flags */
typedef struct { int card, M; mpq_t* G; mpq_t* X; char* kG; char* kX; } slot;

static slot* sl_new(int card, int M)
{
	slot* s = (slot*) calloc(1, sizeof(slot));
	int i;
	s->card = card; s->M = M;
	s->G = (mpq_t*) malloc(card * sizeof(mpq_t));
	s->X = (mpq_t*) malloc((size_t) card * M * sizeof(mpq_t));
	for (i = 0; i < card; i++) mpq_init(s->G[i]);
	for (i = 0; i < card * M; i++) mpq_init(s->X[i]);
	s->kG = (char*) calloc(card, 1); s->kX = (char*) calloc((size_t) card * M, 1);
	return s;
}

static void sl_free(slot* s)
{
	int i;
	for (i = 0; i < s->card; i++) mpq_clear(s->G[i]);
	for (i = 0; i < s->card * s->M; i++) mpq_clear(s->X[i]);
	free(s->G); free(s->X); free(s->kG); free(s->kX); free(s);
}

static void sl_copy(slot* d, const slot* s)
{
	int i;
	for (i = 0; i < s->card; i++) mpq_set(d->G[i], s->G[i]);
	for (i = 0; i < s->card * s->M; i++) mpq_set(d->X[i], s->X[i]);
	memcpy(d->kG, s->kG, s->card); memcpy(d->kX, s->kX, (size_t) s->card * s->M);
}

/* rows: type 0 = CE (k, c), type 1 = PC (s, c) */
typedef struct { int type, a, c; } zrow;

/* factorization for one pattern of unknowns */
typedef struct {
	int ncol, nrow;          /* nrow = rows kept = ncol */
	int* col;                /* column ids: X: k*card+c ; G: M*card + c */
	zrow* row;               /* chosen rows */
	char* det;               /* col[i] is uniquely determined by the system */
	/* block upper triangular form of the square block (rows matched to
	 * columns, diagonal blocks the strongly connected components, in solve
	 * order): block b holds rows brow[bstart[b] .. bstart[b+1]-1], matched
	 * to columns bcol[...] (indices into row[] and col[]) */
	int nblk;
	int* bstart; int* brow; int* bcol;
	mpz_t** biz;             /* block inverse as integers: inv[t][u] = biz[b][t*bs+u] / bden[b][t] */
	mpz_t** bden;
	int* offn; int** offq; mpq_t** offc;   /* per row: coefficients on columns of earlier blocks */
	int maxbs;
	char* sig; size_t siglen;
} zfact;

static void zf_free(zfact* f)
{
	int i;
	if (!f) return;
	int b, q;
	for (b = 0; b < f->nblk; b++) {
		int bs = f->bstart[b + 1] - f->bstart[b];
		if (f->biz && f->biz[b]) { for (q = 0; q < bs * bs; q++) mpz_clear(f->biz[b][q]); free(f->biz[b]); }
		if (f->bden && f->bden[b]) { for (q = 0; q < bs; q++) mpz_clear(f->bden[b][q]); free(f->bden[b]); }
	}
	if (f->offn) for (i = 0; i < f->ncol; i++) { for (q = 0; q < f->offn[i]; q++) mpq_clear(f->offc[i][q]); free(f->offq[i]); free(f->offc[i]); }
	free(f->biz); free(f->bden); free(f->bstart); free(f->brow); free(f->bcol); free(f->offn); free(f->offq); free(f->offc);
	free(f->det); free(f->col); free(f->row); free(f->sig); free(f);
}

/* memo of corner traces: G(n_1..n_{R'-1}, t), t = 0..n_{R'}, of the network
 * restricted to the first R' classes of the order, keyed by population */
static int zc_failfast = 0;   /* stop an order at the first G neither the PC nor a sub-run gives (OIMOM_ZC_NOFAILFAST) */
static long zc_nsub = 0, zc_budget = 0;   /* corner sub-runs allowed per class order (OIMOM_ZC_SUBRUNS) */
typedef struct zmemo { int R; int N[32]; int ok; mpq_t* tr; struct zmemo* next; } zmemo;

/* shifted export (zc_solve_shift): G(N-b) and G^{+k}(N-b) for n shifts b
 * (reduced class indices), where need says so */
typedef struct { int n; const int* b; const char* needG; const char* needX; mpq_t* G; mpq_t* X; } zshift;

typedef struct {
	int R, M;             /* classes of the reduced problem, stations */
	const int* N;
	mpq_t* Z;             /* R */
	mpq_t** L;            /* M x R */
	const int* mi;
	zmemo** memo;         /* shared by the sub-runs of one top-level run */
	int depth;
	int lev;              /* basis level: decrements |c| <= lev (M in comom) */
	const zshift* sh;     /* shifted export, top level only */
} zprob;

static void zm_free(zmemo* m)
{
	while (m) {
		zmemo* n = m->next; int t;
		if (m->tr) { for (t = 0; t <= m->N[m->R - 1]; t++) mpq_clear(m->tr[t]); free(m->tr); }
		free(m); m = n;
	}
}

static int zc_run(const zprob* P, const int* ord, mpq_t g0, mpq_t* g1, mpq_t** gk, mpq_t* trace);

/* G_t(c) for every t of class ord[idx], from a run on the network restricted
 * to classes ord[0..idx] at population (N_E - c, N_r): the values the class-r
 * PC cannot give because they need an orphan.  NULL if that run fails. */
static mpq_t* zc_corner(const zprob* P, const int* ord, int idx, const int* c)
{
	int i, k, R1 = idx + 1, key[32];
	for (i = 0; i < idx; i++) key[i] = P->N[ord[i]] - c[i];
	key[idx] = P->N[ord[idx]];
	zmemo* m;
	for (m = *P->memo; m; m = m->next)
		if (m->R == R1 && !memcmp(m->N, key, R1 * sizeof(int))) return m->ok ? m->tr : NULL;
	m = (zmemo*) calloc(1, sizeof(zmemo));
	m->R = R1; memcpy(m->N, key, R1 * sizeof(int));
	m->tr = (mpq_t*) malloc((key[idx] + 1) * sizeof(mpq_t));
	for (i = 0; i <= key[idx]; i++) mpq_init(m->tr[i]);
	if (P->depth < 64 && zc_nsub < zc_budget) {
		mpq_t Zs[32]; mpq_t** Ls = (mpq_t**) malloc(P->M * sizeof(mpq_t*));
		int o[32];
		for (i = 0; i < R1; i++) { mpq_init(Zs[i]); mpq_set(Zs[i], P->Z[ord[i]]); o[i] = i; }
		for (k = 0; k < P->M; k++) { Ls[k] = (mpq_t*) malloc(R1 * sizeof(mpq_t)); for (i = 0; i < R1; i++) { mpq_init(Ls[k][i]); mpq_set(Ls[k][i], P->L[k][ord[i]]); } }
		zprob S = { R1, P->M, key, Zs, Ls, P->mi, P->memo, P->depth + 1, P->lev, NULL };
		zc_nsub++;
		m->ok = zc_run(&S, o, NULL, NULL, NULL, m->tr) == 0;
		for (i = 0; i < R1; i++) mpq_clear(Zs[i]);
		for (k = 0; k < P->M; k++) { for (i = 0; i < R1; i++) mpq_clear(Ls[k][i]); free(Ls[k]); }
		free(Ls);
	}
	m->next = *P->memo; *P->memo = m;
	return m->ok ? m->tr : NULL;
}

static int zc_check = -1, zc_debug = 0;

/* the coefficient row of r over the columns, as a dense vector indexed by
 * column position (cpos[id], -1 if not a column); returns 0 if the row has
 * a term on a non-column unknown (unusable). */
static int build_row(const zprob* P, const int* ord, const cset* C, const slot* cur, const zrow* w,
                     const int* cpos, mpq_t* out, int ncol)
{
	int M = P->M, card = C->card, nE = C->nE, i, k;
	for (i = 0; i < ncol; i++) mpq_set_ui(out[i], 0, 1);
	int c = w->c;
	if (w->type == 0) {
		k = w->a;
		int id = k * card + c;
		if (cpos[id] < 0) return 0;
		mpq_set_ui(out[cpos[id]], 1, 1);
		if (!cur->kG[c]) { if (cpos[M * card + c] < 0) return 0; mpq_set_si(out[cpos[M * card + c]], -1, 1); }
		for (i = 0; i < nE; i++) {
			if (!mpq_sgn(P->L[k][ord[i]])) continue;
			int d = C->nb[c * nE + i];
			if (C->imp[d]) continue;
			if (cpos[k * card + d] < 0) return 0;
			mpq_sub(out[cpos[k * card + d]], out[cpos[k * card + d]], P->L[k][ord[i]]);
		}
	} else {
		int s = w->a, r = ord[s];
		int d = C->nb[c * nE + s];
		int a = P->N[r] - C->cv[c * nE + s];
		if (!cur->kG[c] && a) { if (cpos[M * card + c] < 0) return 0; mpq_set_si(out[cpos[M * card + c]], a, 1); }
		if (!C->imp[d]) {
			if (mpq_sgn(P->Z[r]) && !cur->kG[d]) { if (cpos[M * card + d] < 0) return 0; mpq_neg(out[cpos[M * card + d]], P->Z[r]); }
			for (k = 0; k < M; k++) {
				if (!mpq_sgn(P->L[k][r])) continue;
				if (cpos[k * card + d] < 0) return 0;
				mpq_t t; mpq_init(t); mpq_set_si(t, P->mi[k], 1); mpq_mul(t, t, P->L[k][r]);
				mpq_sub(out[cpos[k * card + d]], out[cpos[k * card + d]], t); mpq_clear(t);
			}
		}
	}
	return 1;
}

/* right-hand side of row w from the known values */
static void row_rhs(const zprob* P, const int* ord, int r, const cset* C, const slot* cur, const slot* prev,
                    const zrow* w, mpq_t rhs)
{
	int card = C->card, nE = C->nE, c = w->c;
	mpq_t t; mpq_init(t);
	mpq_set_ui(rhs, 0, 1);
	if (w->type == 0) {
		int k = w->a;
		if (mpq_sgn(P->L[k][r])) { mpq_mul(t, P->L[k][r], prev->X[k * card + c]); mpq_add(rhs, rhs, t); }
		if (cur->kG[c]) mpq_add(rhs, rhs, cur->G[c]);
	} else {
		int s = w->a, rs = ord[s];
		int d = C->nb[c * nE + s];
		int a = P->N[rs] - C->cv[c * nE + s];
		if (cur->kG[c] && a) { mpq_set_si(t, a, 1); mpq_mul(t, t, cur->G[c]); mpq_sub(rhs, rhs, t); }
		if (!C->imp[d] && mpq_sgn(P->Z[rs]) && cur->kG[d]) { mpq_mul(t, P->Z[rs], cur->G[d]); mpq_add(rhs, rhs, t); }
	}
	mpq_clear(t);
}

/* Rows usable at this step and the columns they touch. */
static int list_rows(const zprob* P, const int* ord, int r, const cset* C, const slot* cur, const slot* prev,
                     zrow* rows)
{
	int M = P->M, card = C->card, nE = C->nE, n = 0, c, k, i;
	for (c = 0; c < card; c++) {
		if (C->imp[c]) continue;
		for (k = 0; k < M; k++) {
			int ok = 1;
			for (i = 0; i < nE; i++)
				if (mpq_sgn(P->L[k][ord[i]]) && C->nb[c * nE + i] < 0) { ok = 0; break; }
			if (!ok) continue;
			if (mpq_sgn(P->L[k][r]) && !prev->kX[k * card + c]) continue;
			rows[n].type = 0; rows[n].a = k; rows[n].c = c; n++;
		}
		if (C->sum[c] < C->M)
			for (i = 0; i < nE; i++) { rows[n].type = 1; rows[n].a = i; rows[n].c = c; n++; }
	}
	return n;
}

/* maximum matching of rows to columns on the pattern of an n x n matrix
 * (Kuhn's augmenting paths) */
static int bt_aug(int i, const int* ap, const int* aj, int* mcol, int* mrow, char* seen)
{
	int q;
	for (q = ap[i]; q < ap[i + 1]; q++) {
		int j = aj[q];
		if (seen[j]) continue;
		seen[j] = 1;
		if (mrow[j] < 0 || bt_aug(mrow[j], ap, aj, mcol, mrow, seen)) { mcol[i] = j; mrow[j] = i; return 1; }
	}
	return 0;
}

/* Tarjan's strongly connected components on rows: i -> i' when row i has a
 * term on the column matched to i'.  Components come out sinks first, which
 * is the solve order. */
typedef struct { const int* ap; const int* aj; const int* mrow; int* idx; int* low; int* stk; char* on; int ns, cnt; int* comp; int ncomp; int* order; int no; } bt_tarjan;
static void bt_scc(bt_tarjan* T, int i)
{
	int q;
	T->idx[i] = T->low[i] = T->cnt++;
	T->stk[T->ns++] = i; T->on[i] = 1;
	for (q = T->ap[i]; q < T->ap[i + 1]; q++) {
		int w = T->mrow[T->aj[q]];
		if (w == i) continue;
		if (T->idx[w] < 0) { bt_scc(T, w); if (T->low[w] < T->low[i]) T->low[i] = T->low[w]; }
		else if (T->on[w] && T->idx[w] < T->low[i]) T->low[i] = T->idx[w];
	}
	if (T->low[i] == T->idx[i]) {
		int w;
		do { w = T->stk[--T->ns]; T->on[w] = 0; T->comp[w] = T->ncomp; T->order[T->no++] = w; } while (w != i);
		T->ncomp++;
	}
}

/* Build and factor the system for the current pattern: pivot rows and columns,
 * the block triangular form of their square block with the inverse of each
 * diagonal block, and which pivot unknowns are determined. */
static zfact* zf_build(const zprob* P, const int* ord, int r, const cset* C, const slot* cur, const slot* prev,
                       zrow* rows, int nrows, int* cpos)
{
	int M = P->M, card = C->card, ncand = card * (M + 1), i, j, k;
	/* candidate columns: X of possible c, and G of possible c that neither the
	 * class-r PC nor a corner run gave */
	for (i = 0; i < ncand; i++) cpos[i] = -1;
	int nc = 0;
	int* cid = (int*) malloc(ncand * sizeof(int));
	for (i = 0; i < card; i++) {
		if (C->imp[i]) continue;
		for (k = 0; k < M; k++) { cpos[k * card + i] = nc; cid[nc++] = k * card + i; }
		if (!cur->kG[i]) { cpos[M * card + i] = nc; cid[nc++] = M * card + i; }
	}
	/* dense rows over candidates; drop unusable rows, then unused columns */
	mpq_t* A = (mpq_t*) malloc((size_t) nrows * nc * sizeof(mpq_t));
	for (i = 0; i < nrows * nc; i++) mpq_init(A[i]);
	int* rok = (int*) malloc(nrows * sizeof(int));
	char* used = (char*) calloc(nc, 1);
	for (i = 0; i < nrows; i++) {
		rok[i] = build_row(P, ord, C, cur, &rows[i], cpos, A + (size_t) i * nc, nc);
		if (rok[i]) { int any = 0; for (j = 0; j < nc; j++) if (mpq_sgn(A[(size_t) i * nc + j])) { used[j] = 1; any = 1; } if (!any) rok[i] = 0; }
	}
	int nu = 0;
	int* upos = (int*) malloc(nc * sizeof(int));
	for (j = 0; j < nc; j++) upos[j] = used[j] ? nu++ : -1;
	/* forward elimination on a copy: pivot rows and pivot columns */
	int nr = 0;
	int* rid = (int*) malloc(nrows * sizeof(int));
	for (i = 0; i < nrows; i++) if (rok[i]) rid[nr++] = i;
	int* ucid = (int*) malloc((nu ? nu : 1) * sizeof(int));
	for (j = 0; j < nc; j++) if (upos[j] >= 0) ucid[upos[j]] = cid[j];
	mpq_t* B = (mpq_t*) malloc((size_t) (nr ? nr : 1) * (nu ? nu : 1) * sizeof(mpq_t));
	for (i = 0; i < nr; i++) for (j = 0; j < nc; j++) if (upos[j] >= 0) { mpq_init(B[(size_t) i * nu + upos[j]]); mpq_set(B[(size_t) i * nu + upos[j]], A[(size_t) rid[i] * nc + j]); }
	int* rperm = (int*) malloc((nr ? nr : 1) * sizeof(int));
	for (i = 0; i < nr; i++) rperm[i] = i;
	int* pivc = (int*) malloc((nu ? nu : 1) * sizeof(int));
	int* ppos = (int*) malloc((nu ? nu : 1) * sizeof(int));   /* column -> pivot index, -1 free */
	int sing = 0, col, np = 0;
	mpq_t f, t; mpq_init(f); mpq_init(t);
	for (col = 0; col < nu; col++) {
		int p = -1;
		ppos[col] = -1;
		for (i = np; i < nr; i++) if (mpq_sgn(B[(size_t) rperm[i] * nu + col])) { p = i; break; }
		if (p < 0) continue;
		int tmp = rperm[np]; rperm[np] = rperm[p]; rperm[p] = tmp;
		mpq_t* pr = B + (size_t) rperm[np] * nu;
		for (i = np + 1; i < nr; i++) {
			mpq_t* ri = B + (size_t) rperm[i] * nu;
			if (!mpq_sgn(ri[col])) continue;
			mpq_div(f, ri[col], pr[col]);
			for (j = col; j < nu; j++) if (mpq_sgn(pr[j])) { mpq_mul(t, f, pr[j]); mpq_sub(ri[j], ri[j], t); }
		}
		pivc[np] = col; ppos[col] = np; np++;
	}
	/* null space: one vector per free column, by back substitution on the
	 * echelon rows; an unknown is determined iff every null vector is 0 there */
	char* und = (char*) calloc(nu ? nu : 1, 1);
	int nfree = 0;
	mpq_t* x = (mpq_t*) malloc((nu ? nu : 1) * sizeof(mpq_t));
	for (j = 0; j < nu; j++) mpq_init(x[j]);
	for (col = 0; col < nu; col++) {
		if (ppos[col] >= 0) continue;
		nfree++;
		for (j = 0; j < nu; j++) mpq_set_ui(x[j], j == col, 1);
		for (i = np - 1; i >= 0; i--) {
			mpq_t* pr = B + (size_t) rperm[i] * nu;
			int pc = pivc[i];
			mpq_set_ui(f, 0, 1);
			for (j = pc + 1; j < nu; j++) if (mpq_sgn(pr[j]) && mpq_sgn(x[j])) { mpq_mul(t, pr[j], x[j]); mpq_add(f, f, t); }
			mpq_div(x[pc], f, pr[pc]); mpq_neg(x[pc], x[pc]);
		}
		for (j = 0; j < nu; j++) if (mpq_sgn(x[j])) und[j] = 1;
	}
	for (j = 0; j < nu; j++) mpq_clear(x[j]);
	free(x);
	if (zc_debug > 1 && P->depth == 0 && nfree) {
		int nund = 0; for (j = 0; j < nu; j++) nund += und[j];
		fprintf(stderr, "zc:   order");
		for (j = 0; j <= C->nE; j++) fprintf(stderr, " %d", ord[j]);
		fprintf(stderr, ": rank %d of %d, %d unknown(s) undetermined\n", np, nu, nund);
	}
	/* shifted export: a step with every unknown undetermined loses G(N) for
	 * good (no G^{+k} at this step, so no PC value at the next): drop the order */
	if (P->sh && nu > 0) {
		int nund = 0; for (j = 0; j < nu; j++) nund += und[j];
		if (nund == nu) sing = 1;
	}
	for (i = 0; i < nr * nu; i++) mpq_clear(B[i]);
	free(B);
	zfact* F = NULL;
	{
		int n = np;
		F = (zfact*) calloc(1, sizeof(zfact));
		F->ncol = F->nrow = n;
		F->col = (int*) malloc((n ? n : 1) * sizeof(int));
		F->det = (char*) malloc(n ? n : 1);
		F->row = (zrow*) malloc((n ? n : 1) * sizeof(zrow));
		for (i = 0; i < n; i++) { F->col[i] = ucid[pivc[i]]; F->det[i] = !und[pivc[i]]; }
		/* square block: the pivot rows over the pivot columns (free ones at 0),
		 * sparse by rows */
		int nu_save = nu;
		int* ap = (int*) calloc(n + 1, sizeof(int));
		for (i = 0; i < n; i++) {
			int ri = rid[rperm[i]];
			F->row[i] = rows[ri];
			ap[i + 1] = ap[i];
			for (j = 0; j < nc; j++) if (upos[j] >= 0 && ppos[upos[j]] >= 0 && mpq_sgn(A[(size_t) ri * nc + j])) ap[i + 1]++;
		}
		int* aj = (int*) malloc((ap[n] ? ap[n] : 1) * sizeof(int));
		mpq_t* ax = (mpq_t*) malloc((ap[n] ? ap[n] : 1) * sizeof(mpq_t));
		for (i = 0; i < n; i++) {
			int ri = rid[rperm[i]], q = ap[i];
			for (j = 0; j < nc; j++) if (upos[j] >= 0 && ppos[upos[j]] >= 0 && mpq_sgn(A[(size_t) ri * nc + j])) {
				aj[q] = ppos[upos[j]]; mpq_init(ax[q]); mpq_set(ax[q], A[(size_t) ri * nc + j]); q++;
			}
		}
		/* perfect matching (the block is nonsingular), then components */
		int* mcol = (int*) malloc((n ? n : 1) * sizeof(int));
		int* mrow = (int*) malloc((n ? n : 1) * sizeof(int));
		char* seen = (char*) malloc(n ? n : 1);
		for (i = 0; i < n; i++) mcol[i] = mrow[i] = -1;
		for (i = 0; i < n && !sing; i++) { memset(seen, 0, n); if (!bt_aug(i, ap, aj, mcol, mrow, seen)) sing = 1; }
		bt_tarjan T = { ap, aj, mrow, NULL, NULL, NULL, NULL, 0, 0, NULL, 0, NULL, 0 };
		T.idx = (int*) malloc((n ? n : 1) * sizeof(int)); T.low = (int*) malloc((n ? n : 1) * sizeof(int));
		T.stk = (int*) malloc((n ? n : 1) * sizeof(int)); T.on = (char*) calloc(n ? n : 1, 1);
		T.comp = (int*) malloc((n ? n : 1) * sizeof(int)); T.order = (int*) malloc((n ? n : 1) * sizeof(int));
		for (i = 0; i < n; i++) T.idx[i] = -1;
		if (!sing) for (i = 0; i < n; i++) if (T.idx[i] < 0) bt_scc(&T, i);
		F->nblk = sing ? 0 : T.ncomp;
		F->bstart = (int*) calloc(F->nblk + 1, sizeof(int));
		F->brow = (int*) malloc((n ? n : 1) * sizeof(int));
		F->bcol = (int*) malloc((n ? n : 1) * sizeof(int));
		F->biz = (mpz_t**) calloc(F->nblk ? F->nblk : 1, sizeof(mpz_t*));
		F->bden = (mpz_t**) calloc(F->nblk ? F->nblk : 1, sizeof(mpz_t*));
		F->offn = (int*) calloc(n ? n : 1, sizeof(int));
		F->offq = (int**) calloc(n ? n : 1, sizeof(int*));
		F->offc = (mpq_t**) calloc(n ? n : 1, sizeof(mpq_t*));
		F->maxbs = 0;
		if (!sing) {
			/* T.order lists rows component by component, sinks first */
			int b2 = 0;
			for (i = 0; i < n; i++) {
				int w = T.order[i];
				F->brow[i] = w; F->bcol[i] = mcol[w];
				if (i == 0 || T.comp[w] != T.comp[T.order[i - 1]]) F->bstart[b2++] = i;
			}
			F->bstart[F->nblk] = n;
			int* lpos = (int*) malloc((n ? n : 1) * sizeof(int));   /* column -> position in its block, -1 outside */
			for (j = 0; j < n; j++) lpos[j] = -1;
			for (b2 = 0; b2 < F->nblk && !sing; b2++) {
				int s0 = F->bstart[b2], bs = F->bstart[b2 + 1] - s0, u, q;
				if (bs > F->maxbs) F->maxbs = bs;
				for (u = 0; u < bs; u++) lpos[F->bcol[s0 + u]] = u;
				/* off-block coefficients, and the dense diagonal block next to the identity */
				int w2 = 2 * bs;
				mpq_t* S = (mpq_t*) malloc((size_t) bs * w2 * sizeof(mpq_t));
				for (u = 0; u < bs * w2; u++) mpq_init(S[u]);
				for (u = 0; u < bs; u++) {
					int w = F->brow[s0 + u], no = 0;
					for (q = ap[w]; q < ap[w + 1]; q++) if (lpos[aj[q]] < 0) no++;
					F->offn[w] = no;
					F->offq[w] = (int*) malloc((no ? no : 1) * sizeof(int));
					F->offc[w] = (mpq_t*) malloc((no ? no : 1) * sizeof(mpq_t));
					no = 0;
					for (q = ap[w]; q < ap[w + 1]; q++) {
						if (lpos[aj[q]] < 0) { F->offq[w][no] = aj[q]; mpq_init(F->offc[w][no]); mpq_set(F->offc[w][no], ax[q]); no++; }
						else mpq_set(S[(size_t) u * w2 + lpos[aj[q]]], ax[q]);
					}
					mpq_set_ui(S[(size_t) u * w2 + bs + u], 1, 1);
				}
				for (u = 0; u < bs; u++) lpos[F->bcol[s0 + u]] = -1;
				/* Gauss-Jordan on the block */
				for (col = 0; col < bs && !sing; col++) {
					int p = -1;
					for (i = col; i < bs; i++) if (mpq_sgn(S[(size_t) i * w2 + col])) { p = i; break; }
					if (p < 0) { sing = 1; break; }
					if (p != col) for (j = 0; j < w2; j++) mpq_swap(S[(size_t) p * w2 + j], S[(size_t) col * w2 + j]);
					mpq_t* pr = S + (size_t) col * w2;
					mpq_inv(f, pr[col]);
					for (j = 0; j < w2; j++) if (mpq_sgn(pr[j])) mpq_mul(pr[j], pr[j], f);
					for (i = 0; i < bs; i++) {
						if (i == col) continue;
						mpq_t* ri = S + (size_t) i * w2;
						if (!mpq_sgn(ri[col])) continue;
						mpq_set(f, ri[col]);
						for (j = 0; j < w2; j++) if (mpq_sgn(pr[j])) { mpq_mul(t, f, pr[j]); mpq_sub(ri[j], ri[j], t); }
					}
				}
				/* inverse as integer rows over one denominator each */
				F->biz[b2] = (mpz_t*) malloc((size_t) bs * bs * sizeof(mpz_t));
				F->bden[b2] = (mpz_t*) malloc(bs * sizeof(mpz_t));
				for (i = 0; i < bs; i++) {
					mpz_init_set_ui(F->bden[b2][i], 1);
					for (j = 0; j < bs; j++) mpz_lcm(F->bden[b2][i], F->bden[b2][i], mpq_denref(S[(size_t) i * w2 + bs + j]));
					for (j = 0; j < bs; j++) {
						mpz_t* z = &F->biz[b2][(size_t) i * bs + j];
						mpz_init(*z);
						mpz_divexact(*z, F->bden[b2][i], mpq_denref(S[(size_t) i * w2 + bs + j]));
						mpz_mul(*z, *z, mpq_numref(S[(size_t) i * w2 + bs + j]));
					}
				}
				for (u = 0; u < bs * w2; u++) mpq_clear(S[u]);
				free(S);
			}
			if (sing) F->nblk = b2;   /* zf_free clears the blocks built so far */
			free(lpos);
		}
		if (zc_debug > 1 && P->depth == 0 && !sing) fprintf(stderr, "zc:   %d unknown(s) in %d block(s), largest %d\n", n, F->nblk, F->maxbs);
		for (i = 0; i < ap[n]; i++) mpq_clear(ax[i]);
		free(ap); free(aj); free(ax); free(mcol); free(mrow); free(seen);
		free(T.idx); free(T.low); free(T.stk); free(T.on); free(T.comp); free(T.order);
		nu = nu_save;
		if (sing) { zf_free(F); F = NULL; }
	}
	free(ucid); free(pivc); free(ppos); free(und);
	mpq_clear(f); mpq_clear(t);
	for (i = 0; i < nrows * nc; i++) mpq_clear(A[i]);
	free(A); free(rok); free(used); free(upos); free(rid); free(rperm); free(cid);
	(void) prev; (void) r;
	return F;
}

/* signature of the pattern: known flags of G at this step and of X at the
 * previous one (they decide the rows and the columns) */
static char* make_sig(const cset* C, const slot* cur, const slot* prev, size_t* len)
{
	size_t n = C->card + (size_t) C->card * cur->M;
	char* s = (char*) malloc(n);
	memcpy(s, cur->kG, C->card);
	memcpy(s + C->card, prev->kX, (size_t) C->card * cur->M);
	*len = n;
	return s;
}

/* One run with class order ord (a permutation of 0..R-1).  0 on success. */
static int zc_run(const zprob* P, const int* ord, mpq_t g0, mpq_t* g1, mpq_t** gk, mpq_t* trace)
{
	int R = P->R, M = P->M, idx, i, k, t;
	int NE[32];
	slot** hist = NULL; int nhist = 0; cset* Cprev = NULL;
	slot *last = NULL, *last1 = NULL; cset* Clast = NULL;
	int rc = 0, nfact = 0, ncorner = 0;
	clock_t ck0 = clock(), ckf = 0;
	mpq_t v, w; mpq_init(v); mpq_init(w);
	mpz_t db, acc; mpz_init(db); mpz_init(acc);
	for (idx = 0; idx < R && rc == 0; idx++) {
		int r = ord[idx], nE = idx;
		for (i = 0; i < nE; i++) NE[i] = P->N[ord[i]];
		cset* C = cs_new(nE, P->lev, NE);
		int card = C->card, ring = 2;
		slot* sl[2] = { sl_new(card, M), sl_new(card, M) };
		/* t = 0 */
		slot* s0 = sl[0];
		for (i = 0; i < card; i++) {
			if (C->imp[i]) { s0->kG[i] = 1; for (k = 0; k < M; k++) s0->kX[k * card + i] = 1; continue; }
			if (idx == 0) { mpq_set_ui(s0->G[i], 1, 1); s0->kG[i] = 1; for (k = 0; k < M; k++) { mpq_set_ui(s0->X[k * card + i], 1, 1); s0->kX[k * card + i] = 1; } continue; }
			int* c = C->cv + i * nE;
			int j = c[nE - 1];
			int cp[32]; memcpy(cp, c, (nE - 1) * sizeof(int));
			int pi = cs_find(Cprev, cp);
			slot* h = hist[j];
			mpq_set(s0->G[i], h->G[pi]); s0->kG[i] = h->kG[pi];
			for (k = 0; k < M; k++) { mpq_set(s0->X[k * card + i], h->X[k * Cprev->card + pi]); s0->kX[k * card + i] = h->kX[k * Cprev->card + pi]; }
		}
		/* keep the last M+1 steps for the next class (and two for the export) */
		int keep = P->lev + 1;
		slot** nh = (slot**) calloc(keep, sizeof(slot*));
		if (P->N[r] < keep) { nh[P->N[r]] = sl_new(card, M); sl_copy(nh[P->N[r]], s0); }
		int maxrows = card * (M + nE) + 1;
		zrow* rows = (zrow*) malloc(maxrows * sizeof(zrow));
		int* cpos = (int*) malloc((size_t) card * (M + 1) * sizeof(int));
		zfact* F = NULL;
		mpq_t* rhs = NULL; mpq_t* xs = NULL; mpz_t* bz = NULL; int nrhs = 0;
		mpq_t** corner = (mpq_t**) calloc(card, sizeof(mpq_t*));   /* borrowed from the memo */
		int c0 = -1;
		if (trace && idx == R - 1) {
			int z0[32] = {0};
			c0 = cs_find(C, z0);
			if (!s0->kG[c0]) rc = 1; else mpq_set(trace[0], s0->G[c0]);
		}
		for (t = 1; t <= P->N[r] && rc == 0; t++) {
			slot* prev = sl[(t - 1) % ring];
			slot* cur = sl[t % ring];
			/* G_t from the class-r PC */
			for (i = 0; i < card; i++) {
				if (C->imp[i]) { mpq_set_ui(cur->G[i], 0, 1); cur->kG[i] = 1; continue; }
				int ok = 1;
				mpq_set_ui(v, 0, 1);
				if (mpq_sgn(P->Z[r])) { if (!prev->kG[i]) ok = 0; else { mpq_mul(w, P->Z[r], prev->G[i]); mpq_add(v, v, w); } }
				for (k = 0; k < M && ok; k++) {
					if (!mpq_sgn(P->L[k][r])) continue;
					if (!prev->kX[k * card + i]) { ok = 0; break; }
					mpq_set_si(w, P->mi[k], 1); mpq_mul(w, w, P->L[k][r]); mpq_mul(w, w, prev->X[k * card + i]); mpq_add(v, v, w);
				}
				cur->kG[i] = ok;
				if (ok) { mpq_set_si(w, t, 1); mpq_div(cur->G[i], v, w); }
				else {
					if (!corner[i]) corner[i] = zc_corner(P, ord, idx, C->cv + i * nE);
					if (corner[i]) { mpq_set(cur->G[i], corner[i][t]); cur->kG[i] = 1; ncorner++; }
					else if (zc_failfast) { rc = 1; break; }
				}
			}
			for (i = 0; i < card; i++) for (k = 0; k < M; k++) {
				cur->kX[k * card + i] = C->imp[i];
				if (C->imp[i]) mpq_set_ui(cur->X[k * card + i], 0, 1);
			}
			size_t sl_len; char* sig = make_sig(C, cur, prev, &sl_len);
			if (!F || F->siglen != sl_len || memcmp(F->sig, sig, sl_len)) {
				zf_free(F);
				int nrows = list_rows(P, ord, r, C, cur, prev, rows);
				clock_t ck1 = clock();
				F = zf_build(P, ord, r, C, cur, prev, rows, nrows, cpos);
				ckf += clock() - ck1;
				if (!F) { free(sig); rc = 1; break; }
				F->sig = sig; F->siglen = sl_len;
				nfact++;
				if (F->ncol > nrhs) {
					for (i = 0; i < nrhs; i++) { mpq_clear(rhs[i]); mpq_clear(xs[i]); mpz_clear(bz[i]); }
					free(rhs); free(xs); free(bz); nrhs = F->ncol;
					rhs = (mpq_t*) malloc((nrhs ? nrhs : 1) * sizeof(mpq_t));
					bz = (mpz_t*) malloc((nrhs ? nrhs : 1) * sizeof(mpz_t));
					xs = (mpq_t*) malloc((nrhs ? nrhs : 1) * sizeof(mpq_t));
					for (i = 0; i < nrhs; i++) { mpq_init(rhs[i]); mpq_init(xs[i]); mpz_init(bz[i]); }
				}
			} else free(sig);
			int n = F->ncol, j;
			for (i = 0; i < n; i++) row_rhs(P, ord, r, C, cur, prev, &F->row[i], rhs[i]);
			/* block back-substitution: per block, subtract the terms on the
			 * columns of earlier blocks, then apply the block inverse as
			 * integer dot products over a common denominator */
			int bb;
			for (bb = 0; bb < F->nblk; bb++) {
				int s0 = F->bstart[bb], bs = F->bstart[bb + 1] - s0, u, q;
				for (u = 0; u < bs; u++) {
					int ro = F->brow[s0 + u];
					for (q = 0; q < F->offn[ro]; q++) { mpq_mul(w, F->offc[ro][q], xs[F->offq[ro][q]]); mpq_sub(rhs[ro], rhs[ro], w); }
				}
				mpz_set_ui(db, 1);
				for (u = 0; u < bs; u++) mpz_lcm(db, db, mpq_denref(rhs[F->brow[s0 + u]]));
				for (u = 0; u < bs; u++) { mpq_t* rr = &rhs[F->brow[s0 + u]]; mpz_divexact(bz[u], db, mpq_denref(*rr)); mpz_mul(bz[u], bz[u], mpq_numref(*rr)); }
				for (u = 0; u < bs; u++) {
					mpz_set_ui(acc, 0);
					const mpz_t* row = F->biz[bb] + (size_t) u * bs;
					for (j = 0; j < bs; j++) if (mpz_sgn(bz[j]) && mpz_sgn(row[j])) mpz_addmul(acc, row[j], bz[j]);
					mpq_t* xv = &xs[F->bcol[s0 + u]];
					mpz_set(mpq_numref(*xv), acc);
					mpz_mul(mpq_denref(*xv), F->bden[bb][u], db);
					mpq_canonicalize(*xv);
				}
			}
			for (i = 0; i < n; i++) {
				if (!F->det[i]) continue;
				int id = F->col[i];
				if (id >= M * card) { mpq_set(cur->G[id - M * card], xs[i]); cur->kG[id - M * card] = 1; }
				else { mpq_set(cur->X[id], xs[i]); cur->kX[id] = 1; }
			}
			if (zc_check) {
				/* every usable row must hold (the chosen ones do by construction) */
				int nrows = list_rows(P, ord, r, C, cur, prev, rows), q;
				for (q = 0; q < nrows; q++) {
					zrow* wr = &rows[q];
					int c = wr->c, bad = 0, kk;
					mpq_t lhs; mpq_init(lhs);
					if (wr->type == 0) {
						kk = wr->a;
						if (!cur->kX[kk * card + c] || !cur->kG[c]) { mpq_clear(lhs); continue; }
						mpq_sub(lhs, cur->X[kk * card + c], cur->G[c]);
						for (i = 0; i < nE && !bad; i++) {
							if (!mpq_sgn(P->L[kk][ord[i]])) continue;
							int d = C->nb[c * nE + i];
							if (!cur->kX[kk * card + d]) { bad = 1; break; }
							mpq_mul(w, P->L[kk][ord[i]], cur->X[kk * card + d]); mpq_sub(lhs, lhs, w);
						}
						if (!bad && mpq_sgn(P->L[kk][r])) { mpq_mul(w, P->L[kk][r], prev->X[kk * card + c]); mpq_sub(lhs, lhs, w); }
					} else {
						int s = wr->a, rs = ord[s], d = C->nb[c * nE + s];
						if (!cur->kG[c]) { mpq_clear(lhs); continue; }
						mpq_set_si(w, P->N[rs] - C->cv[c * nE + s], 1); mpq_mul(lhs, w, cur->G[c]);
						if (!C->imp[d]) {
							if (mpq_sgn(P->Z[rs])) { if (!cur->kG[d]) bad = 1; else { mpq_mul(w, P->Z[rs], cur->G[d]); mpq_sub(lhs, lhs, w); } }
							for (kk = 0; kk < M && !bad; kk++) {
								if (!mpq_sgn(P->L[kk][rs])) continue;
								if (!cur->kX[kk * card + d]) { bad = 1; break; }
								mpq_set_si(w, P->mi[kk], 1); mpq_mul(w, w, P->L[kk][rs]); mpq_mul(w, w, cur->X[kk * card + d]); mpq_sub(lhs, lhs, w);
							}
						}
					}
					if (!bad && mpq_sgn(lhs)) { fprintf(stderr, "zc: row %d (type %d) violated at class %d step %d\n", q, wr->type, idx, t); rc = 2; }
					mpq_clear(lhs);
				}
			}
			if (P->N[r] - t < keep) { int j2 = P->N[r] - t; nh[j2] = sl_new(card, M); sl_copy(nh[j2], cur); }
			if (c0 >= 0 && rc == 0) { if (!cur->kG[c0]) rc = 1; else mpq_set(trace[t], cur->G[c0]); }
		}
		free(corner);
		zf_free(F);
		for (i = 0; i < nrhs; i++) { mpq_clear(rhs[i]); mpq_clear(xs[i]); mpz_clear(bz[i]); }
		free(rhs); free(xs); free(bz); free(rows); free(cpos);
		if (idx == R - 1 && rc == 0) {
			Clast = C;
			last = nh[0]; nh[0] = NULL;
			if (P->N[r] >= 1) { last1 = nh[1]; nh[1] = NULL; }
		}
		sl_free(sl[0]); sl_free(sl[1]);
		for (i = 0; i < nhist; i++) if (hist[i]) sl_free(hist[i]);
		free(hist);
		if (Cprev) cs_free(Cprev);
		hist = nh; nhist = keep; Cprev = (idx == R - 1 && rc == 0) ? NULL : C;
	}
	if (rc == 0 && g0) {
		/* export: G(N), G(N-e_s), G^{+k}(N-e_s) */
		int nE = R - 1, rl = ord[R - 1], card = Clast->card, z0[32] = {0};
		int c0 = cs_find(Clast, z0);
		if (!last->kG[c0]) { rc = 1; if (zc_debug > 1 && P->depth == 0) fprintf(stderr, "zc:   export: G(N) unknown\n"); }
		else mpq_set(g0, last->G[c0]);
		for (i = 0; i < R && rc == 0; i++) {
			slot* sv; int ci;
			int r = ord[i];
			if (i < nE) { int e[32] = {0}; e[i] = 1; ci = cs_find(Clast, e); sv = last; }
			else { ci = c0; sv = last1; }
			if (!sv || !sv->kG[ci]) { if (zc_debug > 1 && P->depth == 0) fprintf(stderr, "zc:   export: G(N-e_%d) unknown\n", r); rc = 1; break; }
			mpq_set(g1[r], sv->G[ci]);
			for (k = 0; k < M; k++) {
				if (sv->kX[k * card + ci]) mpq_set(gk[k][r], sv->X[k * card + ci]);
				else if (!mpq_sgn(P->L[k][r])) mpq_set_ui(gk[k][r], 0, 1);
				else { if (zc_debug > 1 && P->depth == 0) fprintf(stderr, "zc:   export: G^{+%d}(N-e_%d) unknown\n", k, r); rc = 1; break; }
			}
		}
		(void) rl;
	}
	if (rc == 0 && P->sh) {
		/* shifted export: G(N-b) is G_t(c) of the last class at t = N_r - b_r,
		 * c = b over the earlier classes; that step is in the history while
		 * b_r <= lev */
		const zshift* sh = P->sh;
		int nE = R - 1, rl = ord[R - 1], card = Clast->card, q, c[32];
		for (q = 0; q < sh->n && rc == 0; q++) {
			const int* b = sh->b + (size_t) q * R;
			int sum = 0, j = b[rl];
			for (i = 0; i < nE; i++) { c[i] = b[ord[i]]; sum += c[i]; }
			slot* sv = j == 0 ? last : j == 1 ? last1 : (j <= P->lev ? hist[j] : NULL);
			int ci = sum <= P->lev ? cs_find(Clast, c) : -1;
			if (!sv || ci < 0) { if (zc_debug > 1) fprintf(stderr, "zc:   shift %d outside the basis\n", q); rc = 1; break; }
			if (sh->needG[q]) {
				if (!sv->kG[ci]) { if (zc_debug > 1) fprintf(stderr, "zc:   shift %d: G unknown\n", q); rc = 1; break; }
				mpq_set(sh->G[q], sv->G[ci]);
			}
			for (k = 0; k < M; k++) {
				if (!sh->needX[(size_t) q * M + k]) continue;
				if (!sv->kX[k * card + ci]) { if (zc_debug > 1) fprintf(stderr, "zc:   shift %d: G^{+%d} unknown\n", q, k); rc = 1; break; }
				mpq_set(sh->X[(size_t) q * M + k], sv->X[k * card + ci]);
			}
		}
	}
	if (last) sl_free(last);
	if (last1) sl_free(last1);
	if (Clast) cs_free(Clast);
	for (i = 0; i < nhist; i++) if (hist[i]) sl_free(hist[i]);
	free(hist);
	if (Cprev) cs_free(Cprev);
	mpq_clear(v); mpq_clear(w); mpz_clear(db); mpz_clear(acc);
	if (zc_debug && P->depth == 0) { fprintf(stderr, "zc: order run rc=%d, %d factorization(s), %d corner value(s), %ld sub-run(s), %.2f s (%.2f s factoring)\n", rc, nfact, ncorner, zc_nsub, (double) (clock() - ck0) / CLOCKS_PER_SEC, (double) ckf / CLOCKS_PER_SEC); zc_nsub = 0; }
	return rc;
}

static int next_perm(int* a, int n)
{
	int i = n - 2, j;
	while (i >= 0 && a[i] >= a[i + 1]) i--;
	if (i < 0) return 0;
	j = n - 1;
	while (a[j] <= a[i]) j--;
	int t = a[i]; a[i] = a[j]; a[j] = t;
	for (i++, j = n - 1; i < j; i++, j--) { t = a[i]; a[i] = a[j]; a[j] = t; }
	return 1;
}

/* Public entry: G(N), g1[r] = G(N-e_r), gk[k][r] = G^{+k}(N-e_r) (0 where
 * L_kr = 0 and the value is not in the basis).  Every N_r >= 1, M >= 1,
 * rational Z and L.  0 on success, nonzero if every class order is singular
 * (the caller then uses CoMoM's tiered route).  *tries gets the number of
 * class orders tried. */
int zc_solve(int R, int M, const int* N, mpq_t* Z, mpq_t** L, const int* mi,
             mpq_t g0, mpq_t* g1, mpq_t** gk, int* tries)
{
	int r, k, Rr = 0, map[32], dly[32], nd = 0;
	if (zc_check < 0) { if (getenv("OIMOM_ZC_NOFAILFAST")) zc_failfast = 0; if (getenv("OIMOM_ZC_SUBRUNS")) zc_budget = atol(getenv("OIMOM_ZC_SUBRUNS")); zc_check = getenv("OIMOM_ZCHECK") != NULL; zc_debug = getenv("OIMOM_DEBUG") ? atoi(getenv("OIMOM_DEBUG")) : 0; }
	if (R > 30) return 1;
	for (r = 0; r < R; r++) {
		for (k = 0; k < M; k++) if (mpq_sgn(L[k][r])) break;
		if (k < M) map[Rr++] = r;
		else { if (!mpq_sgn(Z[r])) return 1; dly[nd++] = r; }
	}
	if (Rr == 0) return 1;
	/* reduced problem over the classes some station serves */
	int Nr[32]; mpq_t Zr[32]; mpq_t** Lr = (mpq_t**) malloc(M * sizeof(mpq_t*));
	for (r = 0; r < Rr; r++) { Nr[r] = N[map[r]]; mpq_init(Zr[r]); mpq_set(Zr[r], Z[map[r]]); }
	for (k = 0; k < M; k++) { Lr[k] = (mpq_t*) malloc(Rr * sizeof(mpq_t)); for (r = 0; r < Rr; r++) { mpq_init(Lr[k][r]); mpq_set(Lr[k][r], L[k][map[r]]); } }
	mpq_t h0; mpq_init(h0);
	mpq_t* h1 = (mpq_t*) malloc(Rr * sizeof(mpq_t));
	mpq_t** hk = (mpq_t**) malloc(M * sizeof(mpq_t*));
	for (r = 0; r < Rr; r++) mpq_init(h1[r]);
	for (k = 0; k < M; k++) { hk[k] = (mpq_t*) malloc(Rr * sizeof(mpq_t)); for (r = 0; r < Rr; r++) mpq_init(hk[k][r]); }
	zmemo* memo = NULL;
	zprob P = { Rr, M, Nr, Zr, Lr, mi, &memo, 0, M, NULL };
	/* class orders: the last one that worked on a problem with the same
	 * classes first (derived networks of one model share their zero and
	 * proportionality structure), then every permutation */
	static int last_ord[32], last_R = 0;
	int ord[32], rc = 1, n = 0;
	if (last_R == Rr) {
		n++;
		zc_nsub = 0;
		rc = zc_run(&P, last_ord, h0, h1, hk, NULL);
		zm_free(memo); memo = NULL;
	}
	for (r = 0; r < Rr; r++) ord[r] = r;
	if (rc != 0 && rc != 2) do {
		if (last_R == Rr && !memcmp(ord, last_ord, Rr * sizeof(int))) continue;
		n++;
		zc_nsub = 0;
		rc = zc_run(&P, ord, h0, h1, hk, NULL);
		zm_free(memo); memo = NULL;
		if (rc == 0) { memcpy(last_ord, ord, Rr * sizeof(int)); last_R = Rr; }
		if (rc == 2) break;
	} while (rc != 0 && n < 120 && next_perm(ord, Rr));
	if (tries) *tries = n;
	if (rc == 0) {
		/* delay-only classes: G = prod Z_r^{N_r}/N_r! G_rest */
		mpq_t df, t; mpq_init(df); mpq_init(t);
		mpz_t f; mpz_init(f);
		mpq_set_ui(df, 1, 1);
		int i;
		for (i = 0; i < nd; i++) {
			r = dly[i];
			mpz_pow_ui(mpq_numref(t), mpq_numref(Z[r]), N[r]);
			mpz_pow_ui(mpq_denref(t), mpq_denref(Z[r]), N[r]);
			mpz_fac_ui(f, N[r]); mpz_mul(mpq_denref(t), mpq_denref(t), f);
			mpq_canonicalize(t);
			mpq_mul(df, df, t);
		}
		mpq_mul(g0, h0, df);
		for (r = 0; r < Rr; r++) {
			mpq_mul(g1[map[r]], h1[r], df);
			for (k = 0; k < M; k++) mpq_mul(gk[k][map[r]], hk[k][r], df);
		}
		for (i = 0; i < nd; i++) {
			r = dly[i];
			mpq_set_ui(t, N[r], 1); mpq_mul(g1[r], g0, t); mpq_div(g1[r], g1[r], Z[r]);
			for (k = 0; k < M; k++) mpq_set_ui(gk[k][r], 0, 1);   /* L_kr = 0 */
		}
		mpq_clear(df); mpq_clear(t); mpz_clear(f);
	}
	for (r = 0; r < Rr; r++) { mpq_clear(Zr[r]); mpq_clear(h1[r]); }
	for (k = 0; k < M; k++) { for (r = 0; r < Rr; r++) { mpq_clear(Lr[k][r]); mpq_clear(hk[k][r]); } free(Lr[k]); free(hk[k]); }
	free(Lr); free(h1); free(hk); mpq_clear(h0);
	return rc;
}

/* Shifted export: for q < nb, G(N - b_q) into G[q] where needG[q], and
 * G^{+k}(N - b_q) into X[q*M+k] where needX[q*M+k].  Every b_q <= N.  The
 * basis level is max(M, the largest decrement over the earlier classes), so
 * every shift is in the basis; 0 on success, nonzero if some needed value is
 * undetermined in every class order tried. */
int zc_solve_shift(int R, int M, const int* N, mpq_t* Z, mpq_t** L, const int* mi,
                   int nb, const int* b, const char* needG, const char* needX,
                   mpq_t* G, mpq_t* X, int* tries)
{
	int r, k, q, Rr = 0, map[32], dly[32], nd = 0, lev = M;
	if (zc_check < 0) { if (getenv("OIMOM_ZC_NOFAILFAST")) zc_failfast = 0; if (getenv("OIMOM_ZC_SUBRUNS")) zc_budget = atol(getenv("OIMOM_ZC_SUBRUNS")); zc_check = getenv("OIMOM_ZCHECK") != NULL; zc_debug = getenv("OIMOM_DEBUG") ? atoi(getenv("OIMOM_DEBUG")) : 0; }
	if (R > 30) return 1;
	for (r = 0; r < R; r++) {
		for (k = 0; k < M; k++) if (mpq_sgn(L[k][r])) break;
		if (k < M) map[Rr++] = r;
		else { if (!mpq_sgn(Z[r])) return 1; dly[nd++] = r; }
	}
	if (Rr == 0) return 1;
	for (q = 0; q < nb; q++) {
		int sum = 0;
		for (r = 0; r < Rr; r++) sum += b[(size_t) q * R + map[r]];
		if (sum > lev) lev = sum;
	}
	int Nr[32]; mpq_t Zr[32]; mpq_t** Lr = (mpq_t**) malloc(M * sizeof(mpq_t*));
	for (r = 0; r < Rr; r++) { Nr[r] = N[map[r]]; mpq_init(Zr[r]); mpq_set(Zr[r], Z[map[r]]); }
	for (k = 0; k < M; k++) { Lr[k] = (mpq_t*) malloc(Rr * sizeof(mpq_t)); for (r = 0; r < Rr; r++) { mpq_init(Lr[k][r]); mpq_set(Lr[k][r], L[k][map[r]]); } }
	int* br = (int*) malloc(((size_t) nb * Rr + 1) * sizeof(int));
	for (q = 0; q < nb; q++) for (r = 0; r < Rr; r++) br[(size_t) q * Rr + r] = b[(size_t) q * R + map[r]];
	zshift sh = { nb, br, needG, needX, G, X };
	zmemo* memo = NULL;
	zprob P = { Rr, M, Nr, Zr, Lr, mi, &memo, 0, lev, &sh };
	static int last_ord[32], last_R = 0;
	int ord[32], rc = 1, n = 0;
	if (last_R == Rr) { n++; zc_nsub = 0; rc = zc_run(&P, last_ord, NULL, NULL, NULL, NULL); zm_free(memo); memo = NULL; }
	for (r = 0; r < Rr; r++) ord[r] = r;
	if (rc != 0 && rc != 2) do {
		if (last_R == Rr && !memcmp(ord, last_ord, Rr * sizeof(int))) continue;
		n++;
		zc_nsub = 0;
		rc = zc_run(&P, ord, NULL, NULL, NULL, NULL);
		zm_free(memo); memo = NULL;
		if (rc == 0) { memcpy(last_ord, ord, Rr * sizeof(int)); last_R = Rr; }
		if (rc == 2) break;
	} while (rc != 0 && n < 120 && next_perm(ord, Rr));
	if (tries) *tries = n;
	if (rc == 0 && nd) {
		/* delay-only classes: every value times prod Z_r^{N_r-b_r}/(N_r-b_r)! */
		mpq_t df, t; mpq_init(df); mpq_init(t);
		mpz_t f; mpz_init(f);
		int i;
		for (q = 0; q < nb; q++) {
			mpq_set_ui(df, 1, 1);
			for (i = 0; i < nd; i++) {
				r = dly[i];
				int e = N[r] - b[(size_t) q * R + r];
				mpz_pow_ui(mpq_numref(t), mpq_numref(Z[r]), e);
				mpz_pow_ui(mpq_denref(t), mpq_denref(Z[r]), e);
				mpz_fac_ui(f, e); mpz_mul(mpq_denref(t), mpq_denref(t), f);
				mpq_canonicalize(t);
				mpq_mul(df, df, t);
			}
			if (needG[q]) mpq_mul(G[q], G[q], df);
			for (k = 0; k < M; k++) if (needX[(size_t) q * M + k]) mpq_mul(X[(size_t) q * M + k], X[(size_t) q * M + k], df);
		}
		mpq_clear(df); mpq_clear(t); mpz_clear(f);
	}
	for (r = 0; r < Rr; r++) mpq_clear(Zr[r]);
	for (k = 0; k < M; k++) { for (r = 0; r < Rr; r++) mpq_clear(Lr[k][r]); free(Lr[k]); }
	free(Lr); free(br);
	return rc;
}
