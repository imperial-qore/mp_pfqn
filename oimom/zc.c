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
	mpq_t* inv;              /* ncol x ncol */
	mpz_t* iz;               /* inv as integers: inv[i][j] = iz[i][j] / den[i] */
	mpz_t* den;
	char* sig; size_t siglen;
} zfact;

static void zf_free(zfact* f)
{
	int i;
	if (!f) return;
	for (i = 0; i < f->ncol * f->ncol; i++) { mpq_clear(f->inv[i]); mpz_clear(f->iz[i]); }
	for (i = 0; i < f->ncol; i++) mpz_clear(f->den[i]);
	free(f->inv); free(f->iz); free(f->den); free(f->det); free(f->col); free(f->row); free(f->sig); free(f);
}

/* memo of corner traces: G(n_1..n_{R'-1}, t), t = 0..n_{R'}, of the network
 * restricted to the first R' classes of the order, keyed by population */
static int zc_failfast = 0;   /* stop an order at the first G neither the PC nor a sub-run gives (OIMOM_ZC_NOFAILFAST) */
static long zc_nsub = 0, zc_budget = 0;   /* corner sub-runs allowed per class order (OIMOM_ZC_SUBRUNS) */
typedef struct zmemo { int R; int N[32]; int ok; mpq_t* tr; struct zmemo* next; } zmemo;

typedef struct {
	int R, M;             /* classes of the reduced problem, stations */
	const int* N;
	mpq_t* Z;             /* R */
	mpq_t** L;            /* M x R */
	const int* mi;
	zmemo** memo;         /* shared by the sub-runs of one top-level run */
	int depth;
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
		zprob S = { R1, P->M, key, Zs, Ls, P->mi, P->memo, P->depth + 1 };
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
		if (C->sum[c] < M)
			for (i = 0; i < nE; i++) { rows[n].type = 1; rows[n].a = i; rows[n].c = c; n++; }
	}
	return n;
}

/* Build and factor the system for the current pattern: pivot rows and columns,
 * the inverse of their block, and which pivot unknowns are determined. */
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
		/* square block: the pivot rows over the pivot columns (free ones at 0) */
		mpq_t* S = (mpq_t*) malloc((size_t) (n ? n : 1) * n * 2 * sizeof(mpq_t));
		int w2 = 2 * n;
		for (i = 0; i < n; i++) {
			int ri = rid[rperm[i]];
			F->row[i] = rows[ri];
			for (j = 0; j < w2; j++) mpq_init(S[(size_t) i * w2 + j]);
			for (j = 0; j < nc; j++) if (upos[j] >= 0 && ppos[upos[j]] >= 0) mpq_set(S[(size_t) i * w2 + ppos[upos[j]]], A[(size_t) ri * nc + j]);
			mpq_set_ui(S[(size_t) i * w2 + n + i], 1, 1);
		}
		int nu_save = nu; nu = n;
		for (col = 0; col < nu && !sing; col++) {
			int p = -1;
			for (i = col; i < nu; i++) if (mpq_sgn(S[(size_t) i * w2 + col])) { p = i; break; }
			if (p < 0) { sing = 1; break; }
			if (p != col) for (j = 0; j < w2; j++) mpq_swap(S[(size_t) p * w2 + j], S[(size_t) col * w2 + j]);
			mpq_t* pr = S + (size_t) col * w2;
			mpq_inv(f, pr[col]);
			for (j = 0; j < w2; j++) if (mpq_sgn(pr[j])) mpq_mul(pr[j], pr[j], f);
			for (i = 0; i < nu; i++) {
				if (i == col) continue;
				mpq_t* ri = S + (size_t) i * w2;
				if (!mpq_sgn(ri[col])) continue;
				mpq_set(f, ri[col]);
				for (j = 0; j < w2; j++) if (mpq_sgn(pr[j])) { mpq_mul(t, f, pr[j]); mpq_sub(ri[j], ri[j], t); }
			}
		}
		/* the chosen rows of the row-swapped block keep their order: the
		 * inverse maps the rhs of F->row[i] (i in the original order) */
		F->inv = (mpq_t*) malloc((size_t) (nu ? nu : 1) * nu * sizeof(mpq_t));
		for (i = 0; i < nu; i++) for (j = 0; j < nu; j++) { mpq_init(F->inv[(size_t) i * nu + j]); mpq_set(F->inv[(size_t) i * nu + j], S[(size_t) i * w2 + nu + j]); }
		for (i = 0; i < nu * w2; i++) mpq_clear(S[i]);
		free(S);
		F->iz = (mpz_t*) malloc((size_t) (nu ? nu : 1) * nu * sizeof(mpz_t));
		F->den = (mpz_t*) malloc((nu ? nu : 1) * sizeof(mpz_t));
		for (i = 0; i < nu; i++) {
			mpz_init_set_ui(F->den[i], 1);
			for (j = 0; j < nu; j++) mpz_lcm(F->den[i], F->den[i], mpq_denref(F->inv[(size_t) i * nu + j]));
			for (j = 0; j < nu; j++) {
				mpz_t* z = &F->iz[(size_t) i * nu + j];
				mpz_init(*z);
				mpz_divexact(*z, F->den[i], mpq_denref(F->inv[(size_t) i * nu + j]));
				mpz_mul(*z, *z, mpq_numref(F->inv[(size_t) i * nu + j]));
			}
		}
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
	mpq_t v, w; mpq_init(v); mpq_init(w);
	mpz_t db, acc; mpz_init(db); mpz_init(acc);
	for (idx = 0; idx < R && rc == 0; idx++) {
		int r = ord[idx], nE = idx;
		for (i = 0; i < nE; i++) NE[i] = P->N[ord[i]];
		cset* C = cs_new(nE, M, NE);
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
		int keep = M + 1;
		slot** nh = (slot**) calloc(keep, sizeof(slot*));
		if (P->N[r] < keep) { nh[P->N[r]] = sl_new(card, M); sl_copy(nh[P->N[r]], s0); }
		int maxrows = card * (M + nE) + 1;
		zrow* rows = (zrow*) malloc(maxrows * sizeof(zrow));
		int* cpos = (int*) malloc((size_t) card * (M + 1) * sizeof(int));
		zfact* F = NULL;
		mpq_t* rhs = NULL; mpz_t* bz = NULL; int nrhs = 0;
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
				F = zf_build(P, ord, r, C, cur, prev, rows, nrows, cpos);
				if (!F) { free(sig); rc = 1; break; }
				F->sig = sig; F->siglen = sl_len;
				nfact++;
				if (F->ncol > nrhs) {
					for (i = 0; i < nrhs; i++) { mpq_clear(rhs[i]); mpz_clear(bz[i]); }
					free(rhs); free(bz); nrhs = F->ncol;
					rhs = (mpq_t*) malloc((nrhs ? nrhs : 1) * sizeof(mpq_t));
					bz = (mpz_t*) malloc((nrhs ? nrhs : 1) * sizeof(mpz_t));
					for (i = 0; i < nrhs; i++) { mpq_init(rhs[i]); mpz_init(bz[i]); }
				}
			} else free(sig);
			int n = F->ncol, j;
			for (i = 0; i < n; i++) row_rhs(P, ord, r, C, cur, prev, &F->row[i], rhs[i]);
			/* x = inv * rhs as integer dot products over common denominators */
			mpz_set_ui(db, 1);
			for (j = 0; j < n; j++) mpz_lcm(db, db, mpq_denref(rhs[j]));
			for (j = 0; j < n; j++) { mpz_divexact(bz[j], db, mpq_denref(rhs[j])); mpz_mul(bz[j], bz[j], mpq_numref(rhs[j])); }
			for (i = 0; i < n; i++) {
				mpz_set_ui(acc, 0);
				const mpz_t* row = F->iz + (size_t) i * n;
				for (j = 0; j < n; j++) if (mpz_sgn(bz[j]) && mpz_sgn(row[j])) mpz_addmul(acc, row[j], bz[j]);
				mpz_set(mpq_numref(v), acc);
				mpz_mul(mpq_denref(v), F->den[i], db);
				mpq_canonicalize(v);
				if (!F->det[i]) continue;
				int id = F->col[i];
				if (id >= M * card) { mpq_set(cur->G[id - M * card], v); cur->kG[id - M * card] = 1; }
				else { mpq_set(cur->X[id], v); cur->kX[id] = 1; }
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
		for (i = 0; i < nrhs; i++) { mpq_clear(rhs[i]); mpz_clear(bz[i]); }
		free(rhs); free(bz); free(rows); free(cpos);
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
	if (last) sl_free(last);
	if (last1) sl_free(last1);
	if (Clast) cs_free(Clast);
	for (i = 0; i < nhist; i++) if (hist[i]) sl_free(hist[i]);
	free(hist);
	if (Cprev) cs_free(Cprev);
	mpq_clear(v); mpq_clear(w); mpz_clear(db); mpz_clear(acc);
	if (zc_debug && P->depth == 0) { fprintf(stderr, "zc: order run rc=%d, %d factorization(s), %d corner value(s), %ld sub-run(s)\n", rc, nfact, ncorner, zc_nsub); zc_nsub = 0; }
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
	zprob P = { Rr, M, Nr, Zr, Lr, mi, &memo, 0 };
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
