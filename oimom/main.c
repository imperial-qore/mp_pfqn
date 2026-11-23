/* oimom - exact MoM-family analysis of closed networks with order-independent
 * (OI) stations whose rate depends only on the set of classes present.
 *
 * An OI station with visit weights v_r and support rates mu_S (S the set of
 * classes with at least one job there) has balance function
 *   Phi(n) = (1/mu_supp(n)) sum_{r: n_r>0} v_r Phi(n - e_r),  Phi(0) = 1,
 * and generating function F(z) = sum_n Phi(n) z^n.  Unrolling the recursion
 * over the supports gives, with theta_{S,r} = v_r/mu_S for r in S and 0
 * otherwise,
 *   F(z) = sum over chains {} < S_1 < ... < S_k (S_i = S_{i-1} + {r_i})
 *          of  prod_i (v_{r_i}/mu_{S_i}) z_{r_i} / (1 - theta_{S_i} . z),
 * and 1/(1 - theta . z) is the generating function of a load-independent
 * single-server queue with demands theta.  So G(N) of the network is a
 * weighted sum, over one chain per OI station, of normalizing constants of
 * load-independent networks in which each OI station is replaced by one
 * "virtual queue" per support of its chain, evaluated at N minus one job of
 * each class of the chain's last support.  Each of those is solved exactly
 * by CoMoM (inner.c).  Throughputs and queue lengths follow from the
 * constants each inner solve exports (G(P-1_r), G^{+k}(P-1_r)):
 *   G(N-1_r)       = sum_c w_c G_c(P_c - 1_r)
 *   Q_kr G(N)      = sum_c w_c m_k L_kr G_c^{+k}(P_c - 1_r)            (LI queue k)
 *   Q_jr G(N)      = sum_c w_c ( [r in S_j] G_c(P_c)
 *                               + sum_{T in chain_j} theta_{T,r} G_c^{+T}(P_c - 1_r) )   (OI station j)
 * the first term of the last line being the one job of class r the chain's
 * shift places at the station.  --conv solves the network instead by
 * convolution of the balance functions over the population lattice (the
 * oracle), --validate compares the two exactly. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include "util.h"
#include "singular.h"

int oimom_inner_solve(qnmodel* q, mpq_t* out);
extern int OIMOM_THREADS;
int zc_solve(int R, int M, const int* N, mpq_t* Z, mpq_t** L, const int* mi,
             mpq_t g0, mpq_t* g1, mpq_t** gk, int* tries);
static int oimom_zc = 1;   /* zero-aware CoMoM (zc.c) first; --comom-inner skips it */

static double oi_cputime(void)
{
	struct rusage u; getrusage(RUSAGE_SELF, &u);
	return u.ru_utime.tv_sec + u.ru_stime.tv_sec + 1e-6 * (u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}

/* ---- model ---- */
typedef struct {
	int R, M, K;
	int* N;
	mpq_t* Z;      /* R */
	mpq_t** L;     /* M x R, LI queues */
	int* mi;       /* M */
	mpq_t** v;     /* K x R, OI visit weights */
	mpq_t** mu;    /* K x 2^R, mu[j][S], S a bitmask of classes (bit r = class r+1) */
} oimodel;

static int read_q(FILE* f, mpq_t q)
{
	if (gmp_fscanf(f, "%Qd", q) != 1) return 0;
	mpq_canonicalize(q);
	return 1;
}

static oimodel* oi_read(const char* path)
{
	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "oimom: cannot open %s\n", path); return NULL; }
	oimodel* m = (oimodel*) calloc(1, sizeof(oimodel));
	int r, k, j, S;
	if (fscanf(f, "%d", &m->R) != 1 || m->R < 1 || m->R > 16) goto bad;
	m->N = (int*) calloc(m->R, sizeof(int));
	for (r = 0; r < m->R; r++) if (fscanf(f, "%d", &m->N[r]) != 1 || m->N[r] < 0) goto bad;
	m->Z = (mpq_t*) malloc(m->R * sizeof(mpq_t));
	for (r = 0; r < m->R; r++) { mpq_init(m->Z[r]); if (!read_q(f, m->Z[r])) goto bad; }
	if (fscanf(f, "%d", &m->M) != 1 || m->M < 0) goto bad;
	m->L = (mpq_t**) malloc((m->M ? m->M : 1) * sizeof(mpq_t*));
	m->mi = (int*) calloc(m->M ? m->M : 1, sizeof(int));
	for (k = 0; k < m->M; k++) {
		if (fscanf(f, "%d", &m->mi[k]) != 1 || m->mi[k] < 1) goto bad;
		m->L[k] = (mpq_t*) malloc(m->R * sizeof(mpq_t));
		for (r = 0; r < m->R; r++) { mpq_init(m->L[k][r]); if (!read_q(f, m->L[k][r])) goto bad; }
	}
	char kw[64];
	if (fscanf(f, "%63s", kw) == 1) {
		if (strcmp(kw, "OI") != 0) { fprintf(stderr, "oimom: expected 'OI K' after the station block, got '%s'\n", kw); goto bad2; }
		if (fscanf(f, "%d", &m->K) != 1 || m->K < 0) goto bad;
		int nS = 1 << m->R;
		m->v = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
		m->mu = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
		for (j = 0; j < m->K; j++) {
			m->v[j] = (mpq_t*) malloc(m->R * sizeof(mpq_t));
			for (r = 0; r < m->R; r++) { mpq_init(m->v[j][r]); if (!read_q(f, m->v[j][r])) goto bad; }
			m->mu[j] = (mpq_t*) malloc(nS * sizeof(mpq_t));
			mpq_init(m->mu[j][0]);
			for (S = 1; S < nS; S++) { mpq_init(m->mu[j][S]); if (!read_q(f, m->mu[j][S])) goto bad; }
		}
	}
	fclose(f);
	/* checks: nonnegative data, positive rate on every support the station can reach */
	for (r = 0; r < m->R; r++) if (mpq_sgn(m->Z[r]) < 0) { fprintf(stderr, "oimom: negative think time\n"); return NULL; }
	for (k = 0; k < m->M; k++) for (r = 0; r < m->R; r++) if (mpq_sgn(m->L[k][r]) < 0) { fprintf(stderr, "oimom: negative demand\n"); return NULL; }
	for (j = 0; j < m->K; j++) {
		int vis = 0;
		for (r = 0; r < m->R; r++) {
			if (mpq_sgn(m->v[j][r]) < 0) { fprintf(stderr, "oimom: negative OI visit weight\n"); return NULL; }
			if (mpq_sgn(m->v[j][r]) > 0) vis |= 1 << r;
		}
		for (S = 1; S < (1 << m->R); S++)
			if ((S & ~vis) == 0 && mpq_sgn(m->mu[j][S]) <= 0) {
				fprintf(stderr, "oimom: OI station %d needs a positive rate on support %d (a set of visiting classes)\n", m->M + j + 1, S);
				return NULL;
			}
	}
	return m;
bad:
	fprintf(stderr, "oimom: malformed model file %s\n", path);
bad2:
	fclose(f);
	return NULL;
}

/* ---- virtual-queue decomposition ---- */
typedef struct { int len; int cls[16]; int sup[16]; int S; mpq_t w; } chain;
typedef struct { int n; chain* c; } chainset;

static void chains_dfs(oimodel* m, int j, chain* cur, chainset* cs, int* cap)
{
	int r;
	if (cs->n == *cap) { *cap *= 2; cs->c = (chain*) realloc(cs->c, *cap * sizeof(chain)); }
	chain* c = &cs->c[cs->n++];
	c->len = cur->len; c->S = cur->S;
	memcpy(c->cls, cur->cls, sizeof c->cls); memcpy(c->sup, cur->sup, sizeof c->sup);
	mpq_init(c->w); mpq_set(c->w, cur->w);
	for (r = 0; r < m->R; r++) {
		if (cur->S >> r & 1) continue;
		if (mpq_sgn(m->v[j][r]) == 0 || m->N[r] == 0) continue;
		chain nx; mpq_init(nx.w);
		nx.len = cur->len + 1; nx.S = cur->S | (1 << r);
		memcpy(nx.cls, cur->cls, sizeof nx.cls); memcpy(nx.sup, cur->sup, sizeof nx.sup);
		nx.cls[cur->len] = r; nx.sup[cur->len] = nx.S;
		mpq_div(nx.w, m->v[j][r], m->mu[j][nx.S]);
		mpq_mul(nx.w, nx.w, cur->w);
		chains_dfs(m, j, &nx, cs, cap);
		mpq_clear(nx.w);
	}
}

/* station of a derived load-independent network: LI queue k, or the
 * virtual queue of support T of OI station j */
typedef struct { int virt; int k; int j; int T; } stdesc;

static void st_demand(oimodel* m, const stdesc* s, int r, mpq_t d)
{
	if (!s->virt) { mpq_set(d, m->L[s->k][r]); return; }
	if (!(s->T >> r & 1)) { mpq_set_ui(d, 0, 1); return; }
	mpq_div(d, m->v[s->j][r], m->mu[s->j][s->T]);
}

typedef struct {
	oimodel* m;
	mpz_t* D;              /* per-class integer scale of the derived networks */
	chainset* cs;          /* per OI station */
	stdesc* st; int nst;   /* current derived network */
	int* shift;
	int* pick;             /* chosen chain per OI station */
	mpq_t w;
	/* accumulators */
	mpq_t G; mpq_t* GN1; mpq_t** QLI; mpq_t** QOI;
	long nconf, nsolve, nzc;
	int err;
} acc;

/* Solve the derived network acc->st at population P (assumed served: every
 * class with jobs has a think time or a station with positive demand).
 * Writes g0 = G(P), g1[r] = G(P-1_r) (0 if P_r = 0) and gk[s][r] =
 * G^{+s}(P-1_r) for every station s of acc->st, all unscaled.  0 on
 * success. */
static int net_solve(acc* a, int* P, mpq_t g0, mpq_t* g1, mpq_t** gk)
{
	oimodel* m = a->m;
	int R = m->R, r, s, ra, Ra = 0;
	int act[16];
	for (r = 0; r < R; r++) { mpq_set_ui(g1[r], 0, 1); for (s = 0; s < a->nst; s++) mpq_set_ui(gk[s][r], 0, 1); }
	for (r = 0; r < R; r++) if (P[r] > 0) act[Ra++] = r;
	if (Ra == 0) { mpq_set_ui(g0, 1, 1); return 0; }
	mpq_t d; mpq_init(d);
	int* keep = (int*) malloc((a->nst ? a->nst : 1) * sizeof(int));
	int Mk = 0;
	for (s = 0; s < a->nst; s++) {
		keep[s] = -1;
		for (ra = 0; ra < Ra; ra++) { st_demand(m, &a->st[s], act[ra], d); if (mpq_sgn(d)) break; }
		if (ra < Ra) keep[s] = Mk++;
	}
	if (Mk == 0) {   /* pure delay */
		mpz_t fz; mpz_init(fz);
		mpq_t t; mpq_init(t);
		mpq_set_ui(g0, 1, 1);
		for (ra = 0; ra < Ra; ra++) {
			r = act[ra];
			mpz_pow_ui(mpq_numref(t), mpq_numref(m->Z[r]), P[r]);
			mpz_pow_ui(mpq_denref(t), mpq_denref(m->Z[r]), P[r]);
			mpz_fac_ui(fz, P[r]); mpz_mul(mpq_denref(t), mpq_denref(t), fz);
			mpq_canonicalize(t);
			mpq_mul(g0, g0, t);
		}
		for (ra = 0; ra < Ra; ra++) {   /* Z^(P-1)/(P-1)! = Z^P/P! * P/Z */
			r = act[ra];
			mpq_set_ui(t, P[r], 1); mpq_mul(g1[r], g0, t); mpq_div(g1[r], g1[r], m->Z[r]);
		}
		mpz_clear(fz); mpq_clear(t); mpq_clear(d); free(keep);
		return 0;
	}
	if (oimom_zc) {
		/* zero-aware CoMoM on the rational data: no scaling, no singular fallback */
		int Pa[16], mis[64], tries = 0, rc;
		mpq_t* Zs = sg_mpq_arr(Ra);
		mpq_t** Ls = (mpq_t**) malloc(Mk * sizeof(mpq_t*));
		mpq_t* h1 = sg_mpq_arr(Ra);
		mpq_t** hk = (mpq_t**) malloc(Mk * sizeof(mpq_t*));
		for (s = 0; s < Mk; s++) { Ls[s] = sg_mpq_arr(Ra); hk[s] = sg_mpq_arr(Ra); }
		for (ra = 0; ra < Ra; ra++) { Pa[ra] = P[act[ra]]; mpq_set(Zs[ra], m->Z[act[ra]]); }
		for (s = 0; s < a->nst; s++) {
			if (keep[s] < 0) continue;
			mis[keep[s]] = a->st[s].virt ? 1 : m->mi[a->st[s].k];
			for (ra = 0; ra < Ra; ra++) st_demand(m, &a->st[s], act[ra], Ls[keep[s]][ra]);
		}
		rc = Mk <= 64 ? zc_solve(Ra, Mk, Pa, Zs, Ls, mis, g0, h1, hk, &tries) : 1;
		if (getenv("OIMOM_DEBUG")) {
			fprintf(stderr, "oimom: zc R=%d M=%d: %s after %d class order(s)\n", Ra, Mk, rc ? "singular" : "solved", tries);
			if (rc && getenv("OIMOM_DEBUG")[0] == '2') {
				fprintf(stderr, "  N:"); for (ra = 0; ra < Ra; ra++) fprintf(stderr, " %d", Pa[ra]);
				fprintf(stderr, "\n  Z:"); for (ra = 0; ra < Ra; ra++) gmp_fprintf(stderr, " %Qd", Zs[ra]);
				fprintf(stderr, "\n");
				for (s = 0; s < Mk; s++) { fprintf(stderr, "  %d", mis[s]); for (ra = 0; ra < Ra; ra++) gmp_fprintf(stderr, " %Qd", Ls[s][ra]); fprintf(stderr, "\n"); }
			}
		}
		if (rc == 0) {
			for (ra = 0; ra < Ra; ra++) {
				mpq_set(g1[act[ra]], h1[ra]);
				for (s = 0; s < a->nst; s++) if (keep[s] >= 0) mpq_set(gk[s][act[ra]], hk[keep[s]][ra]);
			}
			a->nsolve++; a->nzc++;
		}
		for (s = 0; s < Mk; s++) { sg_mpq_arr_free(Ls[s], Ra); sg_mpq_arr_free(hk[s], Ra); }
		free(Ls); free(hk); sg_mpq_arr_free(Zs, Ra); sg_mpq_arr_free(h1, Ra);
		if (rc == 0) { mpq_clear(d); free(keep); return 0; }
	}
	qnmodel* q = sg_qn_new(Ra, Mk);
	for (ra = 0; ra < Ra; ra++) {
		r = act[ra];
		q->N[ra] = P[r];
		mpq_set_z(d, a->D[r]); mpq_mul(d, d, m->Z[r]);
		mpz_set(q->Z[ra], mpq_numref(d));          /* integer by the choice of D */
	}
	for (s = 0; s < a->nst; s++) {
		if (keep[s] < 0) continue;
		q->mi[keep[s]] = a->st[s].virt ? 1 : m->mi[a->st[s].k];
		for (ra = 0; ra < Ra; ra++) {
			st_demand(m, &a->st[s], act[ra], d);
			mpq_t t; mpq_init(t); mpq_set_z(t, a->D[act[ra]]); mpq_mul(d, d, t); mpq_clear(t);
			mpz_set(q->L[keep[s]][ra], mpq_numref(d));
		}
	}
	sg_qn_fix_nt(q);
	int nres = SG_NRES(q);
	mpq_t* out = sg_mpq_arr(nres);
	int rc = oimom_inner_solve(q, out);
	a->nsolve++;
	if (rc == 0) {
		/* unscale: a constant at population n carries prod_r D_r^{n_r} */
		mpq_t sc, sc1; mpq_init(sc); mpq_init(sc1);
		mpz_t p; mpz_init(p);
		mpq_set_ui(sc, 1, 1);
		for (ra = 0; ra < Ra; ra++) { mpz_pow_ui(p, a->D[act[ra]], P[act[ra]]); mpz_mul(mpq_numref(sc), mpq_numref(sc), p); }
		mpq_canonicalize(sc);
		mpq_div(g0, out[0], sc);
		for (ra = 0; ra < Ra; ra++) {
			r = act[ra];
			mpq_set_z(sc1, a->D[r]); mpq_div(sc1, sc, sc1);   /* D^{P-1_r} */
			mpq_div(g1[r], out[1+ra], sc1);
			for (s = 0; s < a->nst; s++)
				if (keep[s] >= 0) mpq_div(gk[s][r], out[SG_IDX_GK(Ra, keep[s], ra)], sc1);
		}
		mpq_clear(sc); mpq_clear(sc1); mpz_clear(p);
	} else
		fprintf(stderr, "oimom: inner CoMoM solve failed on a derived network (R=%d, M=%d)\n", Ra, Mk);
	sg_mpq_arr_free(out, nres);
	sg_qn_del(q);
	mpq_clear(d); free(keep);
	return rc;
}

static void process_config(acc* a)
{
	oimodel* m = a->m;
	int R = m->R, r, s, P[16], served[16];
	mpq_t d; mpq_init(d);
	a->nconf++;
	for (r = 0; r < R; r++) P[r] = m->N[r] - a->shift[r];
	/* a class with jobs and nowhere to be makes G vanish */
	int ndeg = 0, rdeg = -1;
	for (r = 0; r < R; r++) {
		if (P[r] == 0) continue;
		served[r] = mpq_sgn(m->Z[r]) != 0;
		for (s = 0; s < a->nst && !served[r]; s++) { st_demand(m, &a->st[s], r, d); if (mpq_sgn(d)) served[r] = 1; }
		if (!served[r]) { ndeg++; rdeg = r; }
	}
	mpq_t g0, t; mpq_init(g0); mpq_init(t);
	mpq_t* g1 = sg_mpq_arr(R);
	mpq_t** gk = (mpq_t**) malloc((a->nst ? a->nst : 1) * sizeof(mpq_t*));
	for (s = 0; s < a->nst; s++) gk[s] = sg_mpq_arr(R);
	if (ndeg == 0) {
		if (net_solve(a, P, g0, g1, gk) != 0) { a->err = 1; goto done; }
		mpq_mul(t, a->w, g0); mpq_add(a->G, a->G, t);
		int j;
		for (j = 0; j < m->K; j++) {
			chain* c = &a->cs[j].c[a->pick[j]];
			for (r = 0; r < R; r++) if (c->S >> r & 1) mpq_add(a->QOI[j][r], a->QOI[j][r], t);
		}
		for (r = 0; r < R; r++) {
			if (P[r] == 0) continue;
			mpq_mul(t, a->w, g1[r]); mpq_add(a->GN1[r], a->GN1[r], t);
			for (s = 0; s < a->nst; s++) {
				st_demand(m, &a->st[s], r, d);
				if (!mpq_sgn(d)) continue;
				mpq_mul(t, a->w, gk[s][r]); mpq_mul(t, t, d);
				if (!a->st[s].virt) {
					mpq_t mm; mpq_init(mm); mpq_set_ui(mm, m->mi[a->st[s].k], 1); mpq_mul(t, t, mm); mpq_clear(mm);
					mpq_add(a->QLI[a->st[s].k][r], a->QLI[a->st[s].k][r], t);
				} else
					mpq_add(a->QOI[a->st[s].j][r], a->QOI[a->st[s].j][r], t);
			}
		}
	} else if (ndeg == 1 && P[rdeg] == 1) {
		/* G_c(P) = 0, but G_c(P - 1_rdeg) has class rdeg empty; its
		 * G^{+k} terms carry a zero demand on rdeg and vanish */
		P[rdeg] = 0;
		if (net_solve(a, P, g0, g1, gk) != 0) { a->err = 1; goto done; }
		mpq_mul(t, a->w, g0); mpq_add(a->GN1[rdeg], a->GN1[rdeg], t);
	}
done:
	for (s = 0; s < a->nst; s++) sg_mpq_arr_free(gk[s], R);
	free(gk); sg_mpq_arr_free(g1, R);
	mpq_clear(g0); mpq_clear(t); mpq_clear(d);
}

static void enum_configs(acc* a, int j)
{
	oimodel* m = a->m;
	if (a->err) return;
	if (j == m->K) { process_config(a); return; }
	int i, r, l;
	for (i = 0; i < a->cs[j].n; i++) {
		chain* c = &a->cs[j].c[i];
		for (r = 0; r < m->R; r++) if ((c->S >> r & 1) && a->shift[r] + 1 > m->N[r]) break;
		if (r < m->R) continue;
		mpq_t wsave; mpq_init(wsave); mpq_set(wsave, a->w);
		mpq_mul(a->w, a->w, c->w);
		for (r = 0; r < m->R; r++) if (c->S >> r & 1) a->shift[r]++;
		int nsave = a->nst;
		for (l = 0; l < c->len; l++) { stdesc* s = &a->st[a->nst++]; s->virt = 1; s->k = -1; s->j = j; s->T = c->sup[l]; }
		a->pick[j] = i;
		enum_configs(a, j + 1);
		a->nst = nsave;
		for (r = 0; r < m->R; r++) if (c->S >> r & 1) a->shift[r]--;
		mpq_set(a->w, wsave); mpq_clear(wsave);
	}
}

/* exact results: G, X[r], Q[s][r] for the M LI queues then the K OI stations */
typedef struct { mpq_t G; mpq_t* X; mpq_t** Q; long nconf, nsolve, nzc; } oires;

static oires* res_new(oimodel* m)
{
	int r, s, S = m->M + m->K;
	oires* o = (oires*) calloc(1, sizeof(oires));
	mpq_init(o->G);
	o->X = sg_mpq_arr(m->R);
	o->Q = (mpq_t**) malloc((S ? S : 1) * sizeof(mpq_t*));
	for (s = 0; s < S; s++) o->Q[s] = sg_mpq_arr(m->R);
	(void) r;
	return o;
}

static oires* oi_solve_mom(oimodel* m)
{
	int R = m->R, j, k, r, s, maxst = m->M + m->K * R;
	acc a; memset(&a, 0, sizeof a);
	a.m = m;
	/* per-class scale making every derived demand and think time integer */
	a.D = (mpz_t*) malloc(R * sizeof(mpz_t));
	mpq_t th; mpq_init(th);
	for (r = 0; r < R; r++) {
		mpz_init_set_ui(a.D[r], 1);
		mpz_lcm(a.D[r], a.D[r], mpq_denref(m->Z[r]));
		for (k = 0; k < m->M; k++) mpz_lcm(a.D[r], a.D[r], mpq_denref(m->L[k][r]));
		for (j = 0; j < m->K; j++) for (s = 1; s < (1 << R); s++) {
			if (!(s >> r & 1) || mpq_sgn(m->mu[j][s]) <= 0) continue;
			mpq_div(th, m->v[j][r], m->mu[j][s]);
			mpz_lcm(a.D[r], a.D[r], mpq_denref(th));
		}
	}
	mpq_clear(th);
	a.cs = (chainset*) calloc(m->K ? m->K : 1, sizeof(chainset));
	for (j = 0; j < m->K; j++) {
		int cap = 16;
		a.cs[j].c = (chain*) malloc(cap * sizeof(chain));
		chain root; memset(&root, 0, sizeof root); mpq_init(root.w); mpq_set_ui(root.w, 1, 1);
		chains_dfs(m, j, &root, &a.cs[j], &cap);
		mpq_clear(root.w);
	}
	a.st = (stdesc*) malloc((maxst ? maxst : 1) * sizeof(stdesc));
	for (k = 0; k < m->M; k++) { a.st[k].virt = 0; a.st[k].k = k; a.st[k].j = -1; a.st[k].T = 0; }
	a.nst = m->M;
	a.shift = (int*) calloc(R, sizeof(int));
	a.pick = (int*) calloc(m->K ? m->K : 1, sizeof(int));
	mpq_init(a.w); mpq_set_ui(a.w, 1, 1);
	mpq_init(a.G);
	a.GN1 = sg_mpq_arr(R);
	a.QLI = (mpq_t**) malloc((m->M ? m->M : 1) * sizeof(mpq_t*));
	for (k = 0; k < m->M; k++) a.QLI[k] = sg_mpq_arr(R);
	a.QOI = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
	for (j = 0; j < m->K; j++) a.QOI[j] = sg_mpq_arr(R);

	enum_configs(&a, 0);
	if (a.err) return NULL;
	if (mpq_sgn(a.G) == 0) { fprintf(stderr, "oimom: G(N) = 0 (a class with jobs has no station it can visit)\n"); return NULL; }

	oires* o = res_new(m);
	mpq_set(o->G, a.G);
	for (r = 0; r < R; r++) mpq_div(o->X[r], a.GN1[r], a.G);
	for (k = 0; k < m->M; k++) for (r = 0; r < R; r++) mpq_div(o->Q[k][r], a.QLI[k][r], a.G);
	for (j = 0; j < m->K; j++) for (r = 0; r < R; r++) mpq_div(o->Q[m->M + j][r], a.QOI[j][r], a.G);
	o->nconf = a.nconf; o->nsolve = a.nsolve; o->nzc = a.nzc;
	return o;
}

/* ---- single basis: one load-independent network for the whole model ----
 * Over a common denominator, F_j(z) = P_j(z) / D_j(z) with D_j = prod over
 * the supports S of (1 - theta_S . z) and P_j a polynomial of degree |D_j|
 * (each chain term is a monomial times the factors of the supports it
 * skips).  Factors of D_j that divide P_j are cancelled.  With G' the
 * normalizing constant of the network of the LI queues plus one queue per
 * remaining support of every OI station, and p = prod_j P_j,
 *   G(N)       = sum_a p_a G'(N-a)
 *   G(N-e_r)   = sum_a p_a G'(N-a-e_r)
 *   Q_kr G(N)  = m_k L_kr sum_a p_a G'^{+k}(N-a-e_r)                           LI queue k
 *   Q_jr G(N)  = sum_a q^{jr}_a G'(N-a) + sum_{S in D_j} theta_{S,r} sum_a p_a G'^{+S}(N-a-e_r)
 * with q^{jr} = (z_r d/dz_r P_j) prod_{i != j} P_i (the derivative of the
 * numerator; the second sum is the derivative of the denominator).  One
 * zero-aware CoMoM run on G' exports every shifted value (zc_solve_shift).
 * Polynomials live on the box n_r <= min(N_r, sum_j |D_j|): coefficients
 * past N_r multiply G' at a negative population. */
int zc_solve_shift(int R, int M, const int* N, mpq_t* Z, mpq_t** L, const int* mi,
                   int nb, const int* b, const char* needG, const char* needX,
                   mpq_t* G, mpq_t* X, int* tries);
static long sb_nq, sb_nshift, sb_tries, sb_deg, sb_cancel;   /* summary of the last single-basis solve */

typedef struct { int R; int B[16]; long st[16]; long size; } sbox;

static void bx_init(sbox* x, int R, const int* B)
{
	int r;
	x->R = R; x->size = 1;
	for (r = 0; r < R; r++) { x->B[r] = B[r]; x->st[r] = x->size; x->size *= B[r] + 1; }
}
static void bx_vec(const sbox* x, long i, int* n) { int r; for (r = 0; r < x->R; r++) n[r] = (int) (i / x->st[r] % (x->B[r] + 1)); }
static long bx_idx(const sbox* x, const int* n) { long i = 0; int r; for (r = 0; r < x->R; r++) { if (n[r] < 0 || n[r] > x->B[r]) return -1; i += n[r] * x->st[r]; } return i; }

/* P *= (1 - th . z), in place (descending: P[i - e_r] is still the old value) */
static void px_linmul(const sbox* x, mpq_t* P, mpq_t* th)
{
	long i; int r, n[16]; mpq_t t; mpq_init(t);
	for (i = x->size - 1; i >= 0; i--) {
		bx_vec(x, i, n);
		for (r = 0; r < x->R; r++)
			if (n[r] && mpq_sgn(th[r]) && mpq_sgn(P[i - x->st[r]])) { mpq_mul(t, th[r], P[i - x->st[r]]); mpq_sub(P[i], P[i], t); }
	}
	mpq_clear(t);
}
/* P /= (1 - th . z) as a power series, in place (ascending) */
static void px_lindiv(const sbox* x, mpq_t* P, mpq_t* th)
{
	long i; int r, n[16]; mpq_t t; mpq_init(t);
	for (i = 0; i < x->size; i++) {
		bx_vec(x, i, n);
		for (r = 0; r < x->R; r++)
			if (n[r] && mpq_sgn(th[r]) && mpq_sgn(P[i - x->st[r]])) { mpq_mul(t, th[r], P[i - x->st[r]]); mpq_add(P[i], P[i], t); }
	}
	mpq_clear(t);
}
/* C = A * B on the box (C distinct from A and B) */
static void px_mul(const sbox* x, mpq_t* C, mpq_t* A, mpq_t* B)
{
	long i, j; int r, n[16], m2[16], s2[16]; mpq_t t; mpq_init(t);
	for (i = 0; i < x->size; i++) mpq_set_ui(C[i], 0, 1);
	for (i = 0; i < x->size; i++) {
		if (!mpq_sgn(A[i])) continue;
		bx_vec(x, i, n);
		for (j = 0; j < x->size; j++) {
			if (!mpq_sgn(B[j])) continue;
			bx_vec(x, j, m2);
			for (r = 0; r < x->R; r++) { s2[r] = n[r] + m2[r]; }
			long k = bx_idx(x, s2);
			if (k < 0) continue;
			mpq_mul(t, A[i], B[j]); mpq_add(C[k], C[k], t);
		}
	}
	mpq_clear(t);
}

static oires* oi_solve_single(oimodel* m)
{
	int R = m->R, Ra = 0, act[16], r, ra, j, k, S, i;
	long ix;
	for (r = 0; r < R; r++) if (m->N[r] > 0) act[Ra++] = r;
	int nU = 1 << Ra;
	/* supports of each OI station over the active classes it visits */
	int* nD = (int*) calloc(m->K ? m->K : 1, sizeof(int));
	char** inD = (char**) malloc((m->K ? m->K : 1) * sizeof(char*));
	mpq_t*** th = (mpq_t***) malloc((m->K ? m->K : 1) * sizeof(mpq_t**));   /* th[j][U][ra] */
	int dtot = 0;
	for (j = 0; j < m->K; j++) {
		int vis = 0;
		for (ra = 0; ra < Ra; ra++) if (mpq_sgn(m->v[j][act[ra]])) vis |= 1 << ra;
		inD[j] = (char*) calloc(nU, 1);
		th[j] = (mpq_t**) calloc(nU, sizeof(mpq_t*));
		for (S = 1; S < nU; S++) {
			if (S & ~vis) continue;
			int so = 0; for (ra = 0; ra < Ra; ra++) if (S >> ra & 1) so |= 1 << act[ra];
			th[j][S] = sg_mpq_arr(Ra);
			for (ra = 0; ra < Ra; ra++) if (S >> ra & 1) mpq_div(th[j][S][ra], m->v[j][act[ra]], m->mu[j][so]);
			inD[j][S] = 1; nD[j]++;
		}
		dtot += nD[j];
	}
	int B[16], Na[16];
	for (ra = 0; ra < Ra; ra++) { Na[ra] = m->N[act[ra]]; B[ra] = Na[ra] < dtot ? Na[ra] : dtot; }
	sbox bx; bx_init(&bx, Ra, B);
	/* P_j = D_j F_j on the box, F_j from the balance recursion */
	mpq_t** Pj = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
	mpq_t t, u; mpq_init(t); mpq_init(u);
	int n[16], n2[16];
	sb_cancel = 0;
	for (j = 0; j < m->K; j++) {
		mpq_t* P = Pj[j] = sg_mpq_arr(bx.size);
		for (ix = 0; ix < bx.size; ix++) {
			bx_vec(&bx, ix, n);
			int sup = 0;
			if (ix == 0) { mpq_set_ui(P[0], 1, 1); continue; }
			for (ra = 0; ra < Ra; ra++) {
				if (!n[ra]) continue;
				sup |= 1 << act[ra];
				mpq_mul(t, m->v[j][act[ra]], P[ix - bx.st[ra]]); mpq_add(P[ix], P[ix], t);
			}
			if (mpq_sgn(P[ix])) mpq_div(P[ix], P[ix], m->mu[j][sup]);
		}
		for (S = 1; S < nU; S++) if (inD[j][S]) px_linmul(&bx, P, th[j][S]);
		/* cancel: (1 - theta_S z) divides P_j iff the quotient has degree < |D_j| */
		mpq_t* Q = sg_mpq_arr(bx.size);
		int changed = 1;
		while (changed) {
			changed = 0;
			for (S = 1; S < nU; S++) {
				if (!inD[j][S]) continue;
				for (ix = 0; ix < bx.size; ix++) mpq_set(Q[ix], P[ix]);
				px_lindiv(&bx, Q, th[j][S]);
				int ok = 1;
				for (ix = 0; ix < bx.size && ok; ix++) {
					if (!mpq_sgn(Q[ix])) continue;
					bx_vec(&bx, ix, n); int tot = 0; for (ra = 0; ra < Ra; ra++) tot += n[ra];
					if (tot >= nD[j]) ok = 0;
				}
				if (!ok) continue;
				for (ix = 0; ix < bx.size; ix++) mpq_set(P[ix], Q[ix]);
				inD[j][S] = 0; nD[j]--; sb_cancel++; changed = 1;
			}
		}
		sg_mpq_arr_free(Q, bx.size);
	}
	/* p = prod P_j, and prod_{i != j} P_i for the numerator derivatives */
	mpq_t* p = sg_mpq_arr(bx.size);
	mpq_t* tmp = sg_mpq_arr(bx.size);
	mpq_t** others = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
	mpq_set_ui(p[0], 1, 1);
	for (j = 0; j < m->K; j++) {
		others[j] = sg_mpq_arr(bx.size); mpq_set_ui(others[j][0], 1, 1);
		for (i = 0; i < m->K; i++) {
			if (i == j) continue;
			px_mul(&bx, tmp, others[j], Pj[i]);
			for (ix = 0; ix < bx.size; ix++) mpq_swap(others[j][ix], tmp[ix]);
		}
	}
	if (m->K) px_mul(&bx, p, others[0], Pj[0]);
	sb_deg = 0;
	for (ix = 0; ix < bx.size; ix++) if (mpq_sgn(p[ix])) { bx_vec(&bx, ix, n); int tot = 0; for (ra = 0; ra < Ra; ra++) tot += n[ra]; if (tot > sb_deg) sb_deg = tot; }
	/* q^{jr} = (z_r d/dz_r P_j) * others_j */
	mpq_t*** q = (mpq_t***) malloc((m->K ? m->K : 1) * sizeof(mpq_t**));
	for (j = 0; j < m->K; j++) {
		q[j] = (mpq_t**) malloc(Ra * sizeof(mpq_t*));
		for (ra = 0; ra < Ra; ra++) {
			for (ix = 0; ix < bx.size; ix++) { bx_vec(&bx, ix, n); mpq_set_ui(t, n[ra], 1); mpq_mul(tmp[ix], Pj[j][ix], t); }
			q[j][ra] = sg_mpq_arr(bx.size);
			px_mul(&bx, q[j][ra], tmp, others[j]);
		}
	}
	/* stations of G': LI queues and the remaining supports, identical rows merged */
	int maxst = m->M + m->K * nU;
	mpq_t** Ls = (mpq_t**) malloc((maxst ? maxst : 1) * sizeof(mpq_t*));
	int* mis = (int*) calloc(maxst ? maxst : 1, sizeof(int));
	int Mk = 0;
	int* liat = (int*) malloc((m->M ? m->M : 1) * sizeof(int));
	int** vqat = (int**) malloc((m->K ? m->K : 1) * sizeof(int*));
	mpq_t* row = sg_mpq_arr(Ra);
	for (k = 0; k < m->M + m->K * nU; k++) {
		int add = 1;
		if (k < m->M) { for (ra = 0; ra < Ra; ra++) mpq_set(row[ra], m->L[k][act[ra]]); add = m->mi[k]; }
		else {
			j = (k - m->M) / nU; S = (k - m->M) % nU;
			if (S == 0) vqat[j] = (int*) malloc(nU * sizeof(int));
			vqat[j][S] = -1;
			if (!S || !inD[j][S]) continue;
			for (ra = 0; ra < Ra; ra++) mpq_set(row[ra], th[j][S][ra]);
		}
		int nz = 0; for (ra = 0; ra < Ra; ra++) if (mpq_sgn(row[ra])) nz = 1;
		int at = -1;
		if (nz) {
			for (i = 0; i < Mk && at < 0; i++) { for (ra = 0; ra < Ra; ra++) if (!mpq_equal(row[ra], Ls[i][ra])) break; if (ra == Ra) at = i; }
			if (at < 0) { at = Mk++; Ls[at] = sg_mpq_arr(Ra); for (ra = 0; ra < Ra; ra++) mpq_set(Ls[at][ra], row[ra]); }
			mis[at] += add;
		}
		if (k < m->M) liat[k] = at; else vqat[j][S] = at;
	}
	sb_nq = Mk;
	/* live classes: served in G' or with a think time; a dead class (no
	 * queue, no think time) keeps G' nonzero only at shifts b_r = N_r */
	int live[16], nl = 0, dead[16];
	for (ra = 0; ra < Ra; ra++) {
		int sv = mpq_sgn(m->Z[act[ra]]) != 0;
		for (i = 0; i < Mk && !sv; i++) if (mpq_sgn(Ls[i][ra])) sv = 1;
		dead[ra] = !sv;
		if (sv) live[nl++] = ra;
	}
	/* requests over the box widened by one */
	int B2[16];
	for (ra = 0; ra < Ra; ra++) B2[ra] = B[ra] + 1 <= Na[ra] ? B[ra] + 1 : Na[ra];
	sbox bx2; bx_init(&bx2, Ra, B2);
	long* rq = (long*) malloc(bx2.size * sizeof(long));
	for (ix = 0; ix < bx2.size; ix++) rq[ix] = -1;
	int nrq = 0, cap = 64;
	int* rb = (int*) malloc((size_t) cap * Ra * sizeof(int));
	char* nG = (char*) malloc(cap); char* nX = (char*) malloc((size_t) cap * (Mk ? Mk : 1));
#define SB_REQ(vec, wantG, kk) do { \
		long i2 = bx_idx(&bx2, vec), w2 = -1; int ok2 = i2 >= 0; \
		for (int d2 = 0; d2 < Ra && ok2; d2++) if (dead[d2] && (vec)[d2] != Na[d2]) ok2 = 0; \
		if (ok2) { \
			if (rq[i2] < 0) { \
				if (nrq == cap) { cap *= 2; rb = (int*) realloc(rb, (size_t) cap * Ra * sizeof(int)); nG = (char*) realloc(nG, cap); nX = (char*) realloc(nX, (size_t) cap * (Mk ? Mk : 1)); } \
				rq[i2] = nrq; memcpy(rb + (size_t) nrq * Ra, vec, Ra * sizeof(int)); nG[nrq] = 0; memset(nX + (size_t) nrq * (Mk ? Mk : 1), 0, Mk ? Mk : 1); nrq++; \
			} \
			w2 = rq[i2]; \
			if (wantG) nG[w2] = 1; \
			if ((kk) >= 0) nX[w2 * (Mk ? Mk : 1) + (kk)] = 1; \
		} } while (0)
	for (ix = 0; ix < bx.size; ix++) {
		int used = mpq_sgn(p[ix]) != 0;
		for (j = 0; j < m->K && !used; j++) for (ra = 0; ra < Ra; ra++) if (mpq_sgn(q[j][ra][ix])) used = 1;
		if (!used) continue;
		bx_vec(&bx, ix, n);
		SB_REQ(n, 1, -1);
		if (!mpq_sgn(p[ix])) continue;
		for (ra = 0; ra < Ra; ra++) {
			if (n[ra] + 1 > Na[ra]) continue;
			memcpy(n2, n, Ra * sizeof(int)); n2[ra]++;
			SB_REQ(n2, 1, -1);
			for (i = 0; i < Mk; i++) if (mpq_sgn(Ls[i][ra])) SB_REQ(n2, 0, i);
		}
	}
#undef SB_REQ
	sb_nshift = nrq;
	mpq_t* gv = sg_mpq_arr(nrq ? nrq : 1);
	mpq_t* xv = sg_mpq_arr((size_t) (nrq ? nrq : 1) * (Mk ? Mk : 1));
	int rc = 0;
	sb_tries = 0;
	if (nl && Mk) {
		/* G' over the live classes */
		int Nl[16], tries = 0, l;
		mpq_t* Zl = sg_mpq_arr(nl);
		mpq_t** Ll = (mpq_t**) malloc(Mk * sizeof(mpq_t*));
		for (l = 0; l < nl; l++) { Nl[l] = Na[live[l]]; mpq_set(Zl[l], m->Z[act[live[l]]]); }
		for (i = 0; i < Mk; i++) { Ll[i] = sg_mpq_arr(nl); for (l = 0; l < nl; l++) mpq_set(Ll[i][l], Ls[i][live[l]]); }
		int* bl = (int*) malloc(((size_t) nrq * nl + 1) * sizeof(int));
		for (i = 0; i < nrq; i++) for (l = 0; l < nl; l++) bl[(size_t) i * nl + l] = rb[(size_t) i * Ra + live[l]];
		rc = zc_solve_shift(nl, Mk, Nl, Zl, Ll, mis, nrq, bl, nG, nX, gv, xv, &tries);
		sb_tries = tries;
		free(bl);
		for (i = 0; i < Mk; i++) sg_mpq_arr_free(Ll[i], nl);
		free(Ll); sg_mpq_arr_free(Zl, nl);
	} else {
		/* no queue: G' is the delay alone over the live classes (no replicas needed) */
		mpz_t f; mpz_init(f);
		for (i = 0; i < nrq; i++) {
			mpq_set_ui(gv[i], 1, 1);
			for (int l = 0; l < nl; l++) {
				int e = Na[live[l]] - rb[(size_t) i * Ra + live[l]];
				mpq_t* z = &m->Z[act[live[l]]];
				mpz_pow_ui(mpq_numref(t), mpq_numref(*z), e); mpz_pow_ui(mpq_denref(t), mpq_denref(*z), e);
				mpz_fac_ui(f, e); mpz_mul(mpq_denref(t), mpq_denref(t), f); mpq_canonicalize(t);
				mpq_mul(gv[i], gv[i], t);
			}
		}
		mpz_clear(f);
	}
	oires* o = NULL;
	if (rc == 0) {
		int MK = Mk ? Mk : 1;
		mpq_t G; mpq_init(G);
		o = res_new(m);
		for (ix = 0; ix < bx.size; ix++) {
			bx_vec(&bx, ix, n);
			long w = bx_idx(&bx2, n), g = w >= 0 ? rq[w] : -1;
			if (g >= 0) {
				if (mpq_sgn(p[ix])) { mpq_mul(t, p[ix], gv[g]); mpq_add(G, G, t); }
				for (j = 0; j < m->K; j++) for (ra = 0; ra < Ra; ra++)
					if (mpq_sgn(q[j][ra][ix])) { mpq_mul(t, q[j][ra][ix], gv[g]); mpq_add(o->Q[m->M + j][act[ra]], o->Q[m->M + j][act[ra]], t); }
			}
			if (!mpq_sgn(p[ix])) continue;
			for (ra = 0; ra < Ra; ra++) {
				if (n[ra] + 1 > Na[ra]) continue;
				memcpy(n2, n, Ra * sizeof(int)); n2[ra]++;
				long w1 = bx_idx(&bx2, n2), g1 = w1 >= 0 ? rq[w1] : -1;
				if (g1 < 0) continue;
				r = act[ra];
				mpq_mul(t, p[ix], gv[g1]); mpq_add(o->X[r], o->X[r], t);
				for (k = 0; k < m->M; k++) {
					if (!mpq_sgn(m->L[k][r]) || liat[k] < 0) continue;
					mpq_mul(t, p[ix], xv[g1 * MK + liat[k]]); mpq_mul(t, t, m->L[k][r]);
					mpq_set_ui(u, m->mi[k], 1); mpq_mul(t, t, u);
					mpq_add(o->Q[k][r], o->Q[k][r], t);
				}
				for (j = 0; j < m->K; j++) for (S = 1; S < nU; S++) {
					if (!inD[j][S] || !(S >> ra & 1) || vqat[j][S] < 0) continue;
					mpq_mul(t, p[ix], xv[g1 * MK + vqat[j][S]]); mpq_mul(t, t, th[j][S][ra]);
					mpq_add(o->Q[m->M + j][r], o->Q[m->M + j][r], t);
				}
			}
		}
		if (mpq_sgn(G) == 0) { fprintf(stderr, "oimom: G(N) = 0 (a class with jobs has no station it can visit)\n"); o = NULL; }
		else {
			mpq_set(o->G, G);
			for (r = 0; r < R; r++) mpq_div(o->X[r], o->X[r], G);
			for (k = 0; k < m->M + m->K; k++) for (r = 0; r < R; r++) mpq_div(o->Q[k][r], o->Q[k][r], G);
			o->nconf = 1; o->nsolve = 1; o->nzc = 1;
		}
		mpq_clear(G);
	}
	/* (buffers are left to process exit, as in oi_solve_mom) */
	mpq_clear(t); mpq_clear(u);
	return o;
}

/* ---- oracle: convolution of the station tables over the population lattice ---- */
typedef struct { int R; int* N; int* rad; long size; } lat;

static void lat_vec(const lat* L, long i, int* n) { int r; for (r = 0; r < L->R; r++) { n[r] = i % (L->N[r] + 1); i /= (L->N[r] + 1); } }
static long lat_idx(const lat* L, const int* n) { long i = 0; int r; for (r = L->R - 1; r >= 0; r--) i = i * (L->N[r] + 1) + n[r]; return i; }

static mpq_t* lat_conv(const lat* L, mpq_t* a, mpq_t* b)
{
	long i, h; int r, n[16], u[16], dlt[16];
	mpq_t* c = sg_mpq_arr(L->size);
	mpq_t t; mpq_init(t);
	for (i = 0; i < L->size; i++) {
		if (!mpq_sgn(a[i])) continue;
		lat_vec(L, i, u);
		for (h = 0; h < L->size; h++) {
			if (!mpq_sgn(b[h])) continue;
			lat_vec(L, h, dlt);
			for (r = 0; r < L->R; r++) { n[r] = u[r] + dlt[r]; if (n[r] > L->N[r]) break; }
			if (r < L->R) continue;
			mpq_mul(t, a[i], b[h]);
			long ix = lat_idx(L, n);
			mpq_add(c[ix], c[ix], t);
		}
	}
	mpq_clear(t);
	return c;
}

static oires* oi_solve_conv(oimodel* m)
{
	int R = m->R, r, s, k, j, n[16];
	long i;
	lat L; L.R = R; L.N = m->N; L.size = 1;
	for (r = 0; r < R; r++) L.size *= (m->N[r] + 1);
	int S = m->M + m->K, F = S + 1;      /* factor 0: delay */
	mpq_t** f = (mpq_t**) malloc(F * sizeof(mpq_t*));
	mpz_t z; mpz_init(z);
	mpq_t t; mpq_init(t);
	for (s = 0; s < F; s++) f[s] = sg_mpq_arr(L.size);
	for (i = 0; i < L.size; i++) {
		lat_vec(&L, i, n);
		int tot = 0; for (r = 0; r < R; r++) tot += n[r];
		/* delay: prod Z^n / n! */
		mpq_set_ui(f[0][i], 1, 1);
		for (r = 0; r < R; r++) {
			mpz_pow_ui(mpq_numref(t), mpq_numref(m->Z[r]), n[r]);
			mpz_pow_ui(mpq_denref(t), mpq_denref(m->Z[r]), n[r]);
			mpz_fac_ui(z, n[r]); mpz_mul(mpq_denref(t), mpq_denref(t), z);
			mpq_canonicalize(t); mpq_mul(f[0][i], f[0][i], t);
		}
		/* LI queue: C(|n|+m-1, m-1) |n|!/prod n_r! prod L^n */
		for (k = 0; k < m->M; k++) {
			mpz_bin_uiui(z, tot + m->mi[k] - 1, m->mi[k] - 1); mpq_set_z(f[1+k][i], z);
			mpz_fac_ui(z, tot); mpq_set_z(t, z); mpq_mul(f[1+k][i], f[1+k][i], t);
			for (r = 0; r < R; r++) {
				mpz_fac_ui(z, n[r]); mpq_set_z(t, z); mpq_div(f[1+k][i], f[1+k][i], t);
				mpz_pow_ui(mpq_numref(t), mpq_numref(m->L[k][r]), n[r]);
				mpz_pow_ui(mpq_denref(t), mpq_denref(m->L[k][r]), n[r]);
				mpq_canonicalize(t); mpq_mul(f[1+k][i], f[1+k][i], t);
			}
		}
		/* OI balance function, in increasing lattice index (n - e_r precedes n) */
		for (j = 0; j < m->K; j++) {
			mpq_t* ph = f[1 + m->M + j];
			if (tot == 0) { mpq_set_ui(ph[i], 1, 1); continue; }
			int sup = 0;
			mpq_set_ui(ph[i], 0, 1);
			for (r = 0; r < R; r++) {
				if (!n[r]) continue;
				sup |= 1 << r;
				n[r]--; mpq_mul(t, m->v[j][r], ph[lat_idx(&L, n)]); n[r]++;
				mpq_add(ph[i], ph[i], t);
			}
			if (mpq_sgn(ph[i])) mpq_div(ph[i], ph[i], m->mu[j][sup]);
		}
	}
	/* prefix and suffix products give every complement G_{-s} */
	mpq_t** pre = (mpq_t**) malloc((F + 1) * sizeof(mpq_t*));
	mpq_t** suf = (mpq_t**) malloc((F + 1) * sizeof(mpq_t*));
	pre[0] = sg_mpq_arr(L.size); mpq_set_ui(pre[0][0], 1, 1);
	suf[F] = sg_mpq_arr(L.size); mpq_set_ui(suf[F][0], 1, 1);
	for (s = 0; s < F; s++) pre[s+1] = lat_conv(&L, pre[s], f[s]);
	for (s = F - 1; s >= 0; s--) suf[s] = lat_conv(&L, suf[s+1], f[s]);
	oires* o = res_new(m);
	long iN = lat_idx(&L, m->N);
	mpq_set(o->G, pre[F][iN]);
	if (mpq_sgn(o->G) == 0) { fprintf(stderr, "oimom: G(N) = 0\n"); return NULL; }
	for (r = 0; r < R; r++) if (m->N[r] > 0) {
		int nm[16]; memcpy(nm, m->N, R * sizeof(int)); nm[r]--;   /* L.N aliases m->N */
		long ix = lat_idx(&L, nm);
		mpq_div(o->X[r], pre[F][ix], o->G);
	}
	int nn[16];
	for (s = 1; s < F; s++) {
		mpq_t* comp = lat_conv(&L, pre[s], suf[s+1]);
		for (i = 0; i < L.size; i++) {
			if (!mpq_sgn(f[s][i])) continue;
			lat_vec(&L, i, n);
			for (r = 0; r < R; r++) nn[r] = m->N[r] - n[r];
			mpq_t* cv = &comp[lat_idx(&L, nn)];
			if (!mpq_sgn(*cv)) continue;
			for (r = 0; r < R; r++) {
				if (!n[r]) continue;
				mpq_mul(t, f[s][i], *cv);
				mpq_t nr; mpq_init(nr); mpq_set_ui(nr, n[r], 1); mpq_mul(t, t, nr); mpq_clear(nr);
				mpq_add(o->Q[s-1][r], o->Q[s-1][r], t);
			}
		}
		sg_mpq_arr_free(comp, L.size);
		for (r = 0; r < R; r++) mpq_div(o->Q[s-1][r], o->Q[s-1][r], o->G);
	}
	mpz_clear(z); mpq_clear(t);
	return o;
}

/* ---- output ---- */
static void put_d(const mpq_t q) { mpf_t f; mpf_init2(f, 128); mpf_set_q(f, q); printf("%.15e", mpf_get_d(f)); mpf_clear(f); }

static int res_equal(oimodel* m, oires* a, oires* b)
{
	int r, s;
	if (!mpq_equal(a->G, b->G)) return 0;
	for (r = 0; r < m->R; r++) if (!mpq_equal(a->X[r], b->X[r])) return 0;
	for (s = 0; s < m->M + m->K; s++) for (r = 0; r < m->R; r++) if (!mpq_equal(a->Q[s][r], b->Q[s][r])) return 0;
	return 1;
}

static void usage(const char* p)
{
	printf("USAGE: %s [-l|-e|-g|-t|-q|-d] [--conv] [--single] [--validate] [--comom-inner] [-j n] model.qnoi\n", p);
	printf("  -l, --log      : log of the normalizing constant\n");
	printf("  -e, --ex       : exact normalizing constant, numerator and denominator\n");
	printf("  -g, --nc       : normalizing constant as double\n");
	printf("  -t, --tput     : throughputs, one per row\n");
	printf("  -q, --qlen     : mean queue lengths, one row per station (LI queues, then OI stations)\n");
	printf("  -d, --exact    : full report in exact rationals\n");
	printf("  --conv         : solve by lattice convolution of the balance functions (oracle)\n");
	printf("  --single       : one network over common denominators (one CoMoM basis) instead of the chain decomposition\n");
	printf("  --validate     : solve both ways and check that every result agrees exactly\n");
	printf("  --comom-inner  : solve every derived network with bin/comom's tiered route instead of\n");
	printf("                   the zero-aware recursion (zc.c) that is tried first\n");
	printf("  -j n           : threads for the per-prime interpolation tier of the inner solves\n");
	printf("Model: a .qn file followed by 'OI K' and K rows 'v_1..v_R mu_1 .. mu_{2^R-1}',\n");
	printf("mu_S the total rate with the set S of classes present (bit r-1 of S = class r).\n");
}

int main(int argc, char** argv)
{
	double t0 = oi_cputime();
	bool lg = false, ex = false, nc = false, tp = false, ql = false, dbg = false, conv = false, val = false, single = false;
	char* file = NULL;
	int i, r, s;
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--log")) lg = true;
		else if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--ex")) ex = true;
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--nc")) nc = true;
		else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tput")) tp = true;
		else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) ql = true;
		else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--exact")) dbg = true;
		else if (!strcmp(argv[i], "--conv")) conv = true;
		else if (!strcmp(argv[i], "--validate")) val = true;
		else if (!strcmp(argv[i], "--single")) single = true;
		else if (!strcmp(argv[i], "--comom-inner")) oimom_zc = 0;
		else if (!strcmp(argv[i], "-j") && i + 1 < argc) OIMOM_THREADS = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
		else if (argv[i][0] == '-') { fprintf(stderr, "oimom: unknown option %s\n", argv[i]); return 1; }
		else file = argv[i];
	}
	if (!file) { usage(argv[0]); return 1; }
	oimodel* m = oi_read(file);
	if (!m) return 1;
	int nt = 0; for (r = 0; r < m->R; r++) nt += m->N[r];
	if (nt == 0) { fprintf(stderr, "oimom: empty network (all populations zero)\n"); return 1; }

	oires* o = NULL;
	if (conv) o = oi_solve_conv(m);
	else if (single) {
		o = oi_solve_single(m);
		if (!o) { fprintf(stderr, "oimom: single basis failed (a shifted value is undetermined), using the chain decomposition\n"); single = false; }
	}
	if (!o && !conv) o = oi_solve_mom(m);
	if (!o) return 1;
	if (val) {
		oires* c = conv ? oi_solve_mom(m) : oi_solve_conv(m);
		if (!c || !res_equal(m, o, c)) { fprintf(stderr, "VALIDATE: FAILED, decomposition and convolution disagree\n"); return 2; }
		fprintf(stderr, "VALIDATE: G, X and Q match the lattice convolution exactly\n");
	}
	int S = m->M + m->K;
	if (lg) printf("%.15e\n", mpq_log(o->G));
	else if (ex) gmp_printf("%Zd\n%Zd\n", mpq_numref(o->G), mpq_denref(o->G));
	else if (nc) { put_d(o->G); printf("\n"); }
	else if (tp) for (r = 0; r < m->R; r++) { put_d(o->X[r]); printf("\n"); }
	else if (ql) for (s = 0; s < S; s++) { for (r = 0; r < m->R; r++) { put_d(o->Q[s][r]); if (r < m->R - 1) printf(" "); } printf("\n"); }
	else {
		printf("Model: R=%d, %d LI queue(s), %d OI station(s)\n", m->R, m->M, m->K);
		if (single) printf("Single basis: one network of %ld queue(s), %ld cancelled factor(s), numerator degree %ld, %ld shifted value(s), %ld class order(s)\n", sb_nq, sb_cancel, sb_deg, sb_nshift, sb_tries);
		else if (!conv) printf("Virtual-queue decomposition: %ld chain configurations, %ld inner CoMoM solves (%ld zero-aware, %ld tiered)\n", o->nconf, o->nsolve, o->nzc, o->nsolve - o->nzc);
		printf("\n========== Performance Metrics ==========\n");
		if (dbg) gmp_printf("G = %Qd\n", o->G);
		else { printf("G = "); put_d(o->G); printf("\nlog(G) = %.15e\n", mpq_log(o->G)); }
		printf("\nX (throughputs):\n");
		for (r = 0; r < m->R; r++) {
			if (dbg) gmp_printf("X[%d] = %Qd\n", r+1, o->X[r]);
			else { printf("X[%d] = ", r+1); put_d(o->X[r]); printf("\n"); }
		}
		printf("\nQ (mean queue lengths):\n");
		mpq_t tot; mpq_init(tot);
		for (s = 0; s < S; s++) {
			printf("Q[%d]%s =", s+1, s >= m->M ? " (OI)" : "");
			mpq_set_ui(tot, 0, 1);
			for (r = 0; r < m->R; r++) {
				mpq_add(tot, tot, o->Q[s][r]);
				if (dbg) gmp_printf("\t%Qd", o->Q[s][r]); else { printf("\t"); put_d(o->Q[s][r]); }
			}
			if (dbg) gmp_printf("\t(total: %Qd)\n", tot); else { printf("\t(total: "); put_d(tot); printf(")\n"); }
		}
		mpq_clear(tot);
		printf("=========================================\n");
		printf("Elapsed time (%s): %.6f s\n", conv ? "OI convolution" : "OIMoM", oi_cputime() - t0);
	}
	return 0;
}
