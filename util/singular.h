#ifndef SINGULAR_H
#define SINGULAR_H

/* Exact fallback for models on which a MoM-type linear system is singular.
 *
 * The solver supplies a callback that returns the exported constants
 *   out[0]               = G(N)
 *   out[1+r]             = G(N-1_r)                  r = 0..R-1
 *   out[SG_IDX_GK(R,k,r)] = G^{+k}(N-1_r)            k = 0..M-1
 * as exact rationals (SG_NRES(q) values), or reports that its system is
 * singular.  sg_solve_tiered() then tries, cheapest first, transformations
 * of the model that leave these constants recoverable exactly:
 *   1. other class orders (up to SG_MAX_CLASS_ORDERS);
 *   2. folding stations with identical demand rows into one station of
 *      summed multiplicity;
 *   3. aggregating classes with proportional demands and think times;
 *   4. the exact limit of a perturbation, by polynomial interpolation.
 * Tiers 1-3 are applied again to each reduced model.  Tier 4 defaults to
 * interpolation over Q from exact solves of perturbed models; a solver
 * may replace it (the multi-modular solvers interpolate per prime).
 *
 * Every tier is exact: none of them approximates G. */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

#define SG_NRES(q) (1 + (q)->R + (q)->M * (q)->R)
#define SG_IDX_GK(R,k,r) (1 + (R) + (k)*(R) + (r))
#define SG_MAX_CLASS_ORDERS 24

/* solve: 0 = solved, out[] filled; 1 = singular system; 2 = hard error
 * (already reported).  interp has the same contract for tier 4; NULL
 * selects sg_interp_exact. */
typedef int (*sg_solve_fn)(qnmodel* q, mpq_t* out, void* ctx);
typedef struct {
	sg_solve_fn solve;
	sg_solve_fn interp;
	void* ctx;
	const char* interp_note;   /* appended to the tier-4 route description */
} sg_ops;

/* ---- derived models ---- */
static inline qnmodel* sg_qn_new(int R, int M)
{
	int r, k;
	qnmodel* q = (qnmodel*) calloc(1, sizeof(qnmodel));
	q->R = R; q->M = M;
	q->N  = (int*) calloc(R, sizeof(int));
	q->mi = (int*) calloc(M, sizeof(int));
	q->Z  = (mpz_t*) malloc(R*sizeof(mpz_t));
	for (r=0;r<R;r++) mpz_init(q->Z[r]);
	q->L  = (mpz_t**) malloc(M*sizeof(mpz_t*));
	for (k=0;k<M;k++) { q->L[k] = (mpz_t*) malloc(R*sizeof(mpz_t)); for (r=0;r<R;r++) mpz_init(q->L[k][r]); }
	return q;
}

static inline void sg_qn_del(qnmodel* q)
{
	int r, k;
	if (!q) return;
	for (k=0;k<q->M;k++) { for (r=0;r<q->R;r++) mpz_clear(q->L[k][r]); free(q->L[k]); }
	for (r=0;r<q->R;r++) mpz_clear(q->Z[r]);
	free(q->L); free(q->Z); free(q->N); free(q->mi); free(q);
}

static inline void sg_qn_fix_nt(qnmodel* q) { int r; q->Nt = 0; for (r=0;r<q->R;r++) q->Nt += q->N[r]; }

static inline mpq_t* sg_mpq_arr(int n) { int i; mpq_t* a = (mpq_t*) malloc(n*sizeof(mpq_t)); for (i=0;i<n;i++) mpq_init(a[i]); return a; }
static inline void sg_mpq_arr_free(mpq_t* a, int n) { int i; if (!a) return; for (i=0;i<n;i++) mpq_clear(a[i]); free(a); }

/* lexicographic next permutation; 0 when the last one was reached */
static inline int sg_next_perm(int* a, int n)
{
	int i = n-2, j, t;
	while (i >= 0 && a[i] >= a[i+1]) i--;
	if (i < 0) return 0;
	j = n-1; while (a[j] <= a[i]) j--;
	t = a[i]; a[i] = a[j]; a[j] = t;
	for (i++, j=n-1; i<j; i++, j--) { t = a[i]; a[i] = a[j]; a[j] = t; }
	return 1;
}

/* Tier 1: the model as given (unless skip_direct, meaning the caller has
 * already seen it fail), then other class orders. */
static inline int sg_solve_order(qnmodel* q, mpq_t* out, int* reordered, const sg_ops* ops, int skip_direct)
{
	int R = q->R, M = q->M, nres = SG_NRES(q), i, j, k, tries = 1, rc = 1;
	*reordered = 0;
	if (!skip_direct) {
		rc = ops->solve(q, out, ops->ctx);
		if (rc != 1) return rc;
	}
	if (R < 2) return rc;
	int* perm = (int*) malloc(R*sizeof(int));
	for (i=0;i<R;i++) perm[i] = i;
	while (rc == 1 && tries < SG_MAX_CLASS_ORDERS && sg_next_perm(perm, R))
	{
		tries++;
		qnmodel* p = sg_qn_new(R, M);
		for (j=0;j<R;j++) {               /* position j holds class perm[j] */
			p->N[j] = q->N[perm[j]];
			mpz_set(p->Z[j], q->Z[perm[j]]);
			for (k=0;k<M;k++) mpz_set(p->L[k][j], q->L[k][perm[j]]);
		}
		for (k=0;k<M;k++) p->mi[k] = q->mi[k];
		sg_qn_fix_nt(p);
		mpq_t* o = sg_mpq_arr(nres);
		rc = ops->solve(p, o, ops->ctx);
		if (rc == 0) {
			mpq_set(out[0], o[0]);
			for (j=0;j<R;j++) {
				mpq_set(out[1+perm[j]], o[1+j]);
				for (k=0;k<M;k++) mpq_set(out[SG_IDX_GK(R,k,perm[j])], o[SG_IDX_GK(R,k,j)]);
			}
			*reordered = 1;
		}
		sg_mpq_arr_free(o, nres);
		sg_qn_del(p);
	}
	free(perm);
	return rc;
}

/* Fold stations with identical demand rows.  A station of multiplicity m
 * is m identical queues, so the folded station has the summed
 * multiplicity; G and G(N-1_r) are unchanged, and an extra replica of any
 * member station is an extra replica of the folded one, so
 * G^{+k}(N-1_r) = G'^{+g(k)}(N-1_r).  Returns NULL if nothing folds. */
static inline qnmodel* sg_merge_stations(qnmodel* q, int* grp)
{
	int M = q->M, R = q->R, k, h, r, nm = 0;
	int* rep = (int*) malloc(M*sizeof(int));
	for (k=0;k<M;k++)
	{
		grp[k] = -1;
		for (h=0;h<nm && grp[k]<0;h++)
		{
			for (r=0;r<R;r++) if (mpz_cmp(q->L[k][r], q->L[rep[h]][r]) != 0) break;
			if (r == R) grp[k] = h;
		}
		if (grp[k] < 0) { rep[nm] = k; grp[k] = nm++; }
	}
	if (nm == M) { free(rep); return NULL; }
	qnmodel* m = sg_qn_new(R, nm);
	for (r=0;r<R;r++) { m->N[r] = q->N[r]; mpz_set(m->Z[r], q->Z[r]); }
	for (h=0;h<nm;h++) for (r=0;r<R;r++) mpz_set(m->L[h][r], q->L[rep[h]][r]);
	for (k=0;k<M;k++) m->mi[grp[k]] += q->mi[k];
	sg_qn_fix_nt(m);
	free(rep);
	return m;
}

/* c with L_kr = c L_ks at every station and Z_r = c Z_s; 0 if none. */
static inline int sg_class_scale(qnmodel* q, int s, int r, mpq_t c)
{
	int k, have = 0;
	mpq_t t; mpq_init(t);
	for (k=0;k<=q->M;k++)
	{
		mpz_srcptr a = (k < q->M) ? q->L[k][s] : q->Z[s];
		mpz_srcptr b = (k < q->M) ? q->L[k][r] : q->Z[r];
		if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) continue;
		if (mpz_sgn(a) == 0 || mpz_sgn(b) == 0) { have = 0; break; }
		mpz_set(mpq_numref(t), b); mpz_set(mpq_denref(t), a); mpq_canonicalize(t);
		if (!have) { mpq_set(c, t); have = 1; }
		else if (!mpq_equal(c, t)) { have = 0; break; }
	}
	mpq_clear(t);
	return have;
}

/* Aggregate classes with proportional demands and think times.  With the
 * 1/n! convention the generating function of G,
 *   exp(sum_r Z_r x_r) prod_k (1 - sum_r L_kr x_r)^(-m_k),
 * depends on such a group only through x_s + sum_i c_i x_i, so
 *   G(N) = G'(n) F(N),  F(N) = prod_groups (sum N_i)!/prod N_i! prod c_i^{N_i},
 * and likewise G(N-1_i) = G'(n-1_g) F(N-1_i) and G^{+k}(N-1_i) =
 * G'^{+k}(n-1_g) F(N-1_i), where ' is the model with each group replaced
 * by its first class at the summed population.  Returns NULL if nothing
 * aggregates; otherwise grp[i] is the merged class of i and c[i] its
 * scale factor. */
static inline qnmodel* sg_merge_classes(qnmodel* q, int* grp, mpq_t* c)
{
	int R = q->R, M = q->M, s, r, k, ng = 0;
	int* rep = (int*) malloc(R*sizeof(int));
	for (r=0;r<R;r++) grp[r] = -1;
	for (s=0;s<R;s++)
	{
		if (grp[s] >= 0) continue;
		rep[ng] = s; grp[s] = ng; mpq_set_ui(c[s], 1, 1);
		for (r=s+1;r<R;r++)
			if (grp[r] < 0 && sg_class_scale(q, s, r, c[r])) grp[r] = ng;
		ng++;
	}
	if (ng == R) { free(rep); return NULL; }
	qnmodel* m = sg_qn_new(ng, M);
	for (s=0;s<ng;s++) {
		mpz_set(m->Z[s], q->Z[rep[s]]);
		for (k=0;k<M;k++) mpz_set(m->L[k][s], q->L[k][rep[s]]);
	}
	for (r=0;r<R;r++) m->N[grp[r]] += q->N[r];
	for (k=0;k<M;k++) m->mi[k] = q->mi[k];
	sg_qn_fix_nt(m);
	free(rep);
	return m;
}

/* F(Nv) of sg_merge_classes, for the group structure grp/c over R classes */
static inline void sg_class_factor(mpq_t F, int* Nv, int R, int ng, int* grp, mpq_t* c)
{
	int g, i;
	mpz_t f; mpz_init(f);
	mpq_t t; mpq_init(t);
	mpq_set_ui(F, 1, 1);
	for (g=0; g<ng; g++)
	{
		unsigned long n = 0;
		for (i=0;i<R;i++) if (grp[i] == g) n += Nv[i];
		mpz_fac_ui(f, n); mpq_set_z(t, f); mpq_mul(F, F, t);
		for (i=0;i<R;i++) if (grp[i] == g)
		{
			mpz_fac_ui(f, Nv[i]); mpq_set_z(t, f); mpq_div(F, F, t);
			mpz_pow_ui(mpq_numref(t), mpq_numref(c[i]), Nv[i]);
			mpz_pow_ui(mpq_denref(t), mpq_denref(c[i]), Nv[i]);
			mpq_canonicalize(t);
			mpq_mul(F, F, t);
		}
	}
	mpz_clear(f); mpq_clear(t);
}

/* Tiers 1-3.  tier receives a short description of the route taken. */
static inline int sg_solve_reduce(qnmodel* q, mpq_t* out, char* tier, size_t tsz, const sg_ops* ops, int skip_direct)
{
	int R = q->R, M = q->M, rc, reordered, i, k;
	rc = sg_solve_order(q, out, &reordered, ops, skip_direct);
	if (rc == 0) { if (reordered) snprintf(tier, tsz, "class order"); return 0; }
	if (rc != 1) return rc;

	/* identical stations */
	{
		int* grp = (int*) malloc(M*sizeof(int));
		qnmodel* ms = sg_merge_stations(q, grp);
		if (ms)
		{
			int nres2 = SG_NRES(ms);
			mpq_t* o = sg_mpq_arr(nres2);
			char sub[96] = "";
			rc = sg_solve_reduce(ms, o, sub, sizeof sub, ops, 0);
			if (rc == 0)
			{
				mpq_set(out[0], o[0]);
				for (i=0;i<R;i++) mpq_set(out[1+i], o[1+i]);
				for (k=0;k<M;k++) for (i=0;i<R;i++) mpq_set(out[SG_IDX_GK(R,k,i)], o[SG_IDX_GK(R,grp[k],i)]);
				snprintf(tier, tsz, "merged stations (M=%d->%d)%s%s", M, ms->M, sub[0] ? ", " : "", sub);
			}
			sg_mpq_arr_free(o, nres2);
			sg_qn_del(ms);
			if (rc != 1) { free(grp); return rc; }
		}
		free(grp);
	}

	/* proportional classes */
	{
		int* grp = (int*) malloc(R*sizeof(int));
		mpq_t* c = sg_mpq_arr(R);
		qnmodel* mc = sg_merge_classes(q, grp, c);
		if (mc)
		{
			int ng = mc->R, nres2 = SG_NRES(mc);
			mpq_t* o = sg_mpq_arr(nres2);
			char sub[96] = "";
			rc = sg_solve_reduce(mc, o, sub, sizeof sub, ops, 0);
			if (rc == 0)
			{
				int* Nv = (int*) malloc(R*sizeof(int));
				mpq_t F; mpq_init(F);
				for (i=0;i<R;i++) Nv[i] = q->N[i];
				sg_class_factor(F, Nv, R, ng, grp, c);
				mpq_mul(out[0], o[0], F);
				for (i=0;i<R;i++)
				{
					if (q->N[i] == 0) {
						mpq_set_ui(out[1+i], 0, 1);
						for (k=0;k<M;k++) mpq_set_ui(out[SG_IDX_GK(R,k,i)], 0, 1);
						continue;
					}
					Nv[i]--;
					sg_class_factor(F, Nv, R, ng, grp, c);
					Nv[i]++;
					mpq_mul(out[1+i], o[1+grp[i]], F);
					for (k=0;k<M;k++) mpq_mul(out[SG_IDX_GK(R,k,i)], o[SG_IDX_GK(ng,k,grp[i])], F);
				}
				mpq_clear(F); free(Nv);
				snprintf(tier, tsz, "merged classes (R=%d->%d)%s%s", R, ng, sub[0] ? ", " : "", sub);
			}
			sg_mpq_arr_free(o, nres2);
			sg_qn_del(mc);
		}
		sg_mpq_arr_free(c, R); free(grp);
	}
	return rc;
}

/* ---- tier 4: exact limit of a perturbation ----
 * Every exported constant is a homogeneous polynomial in (L, Z) jointly,
 * of degree at most Ntot, so along (L, Z) + t (DL, DZ) it is a polynomial
 * g(t) of degree <= Ntot whose value at t=0 is the constant of the
 * original model.  Both L and Z must move: with Z=0 the system can stay
 * singular for every demand matrix.  For a generic direction the system
 * is singular at finitely many t only, so the perturbed model is solved
 * at t=1,2,..., skipping the t that fail, and the Lagrange interpolant
 * through Ntot+1 points is evaluated at t=0.  The direction is drawn from
 * 1..1e6: small integers create coincidences (equal demands) that are
 * themselves MoM degeneracies. */

static inline int sg_interp_degree(const qnmodel* q) { int r, D = 0; for (r=0;r<q->R;r++) D += q->N[r]; return D; }
static inline int sg_interp_points(const qnmodel* q) { return 4*(sg_interp_degree(q)+1) + 8; }

/* The perturbed models at t = 1..T (array index t-1). */
static inline qnmodel** sg_perturb_build(const qnmodel* q, int T)
{
	int R = q->R, M = q->M, r, k, t;
	unsigned long long seed = 0x9E3779B97F4A7C15ULL;
	mpz_t *dz = (mpz_t*) malloc(R*sizeof(mpz_t)), **dl = (mpz_t**) malloc(M*sizeof(mpz_t*));
#define SG_NEXTRAND() (seed = seed*6364136223846793005ULL + 1442695040888963407ULL, 1 + (unsigned long)((seed >> 33) % 1000000))
	for (k=0;k<M;k++) { dl[k] = (mpz_t*) malloc(R*sizeof(mpz_t)); for (r=0;r<R;r++) mpz_init_set_ui(dl[k][r], SG_NEXTRAND()); }
	for (r=0;r<R;r++) mpz_init_set_ui(dz[r], SG_NEXTRAND());
#undef SG_NEXTRAND
	qnmodel** qt = (qnmodel**) malloc(T*sizeof(qnmodel*));
	for (t=1;t<=T;t++)
	{
		qnmodel* m = sg_qn_new(R, M);
		for (r=0;r<R;r++) { m->N[r] = q->N[r]; mpz_mul_ui(m->Z[r], dz[r], t); mpz_add(m->Z[r], m->Z[r], q->Z[r]); }
		for (k=0;k<M;k++) { m->mi[k] = q->mi[k];
			for (r=0;r<R;r++) { mpz_mul_ui(m->L[k][r], dl[k][r], t); mpz_add(m->L[k][r], m->L[k][r], q->L[k][r]); } }
		sg_qn_fix_nt(m);
		qt[t-1] = m;
	}
	for (k=0;k<M;k++) { for (r=0;r<R;r++) mpz_clear(dl[k][r]); free(dl[k]); }
	for (r=0;r<R;r++) mpz_clear(dz[r]);
	free(dl); free(dz);
	return qt;
}

static inline void sg_perturb_free(qnmodel** qt, int T) { int t; for (t=0;t<T;t++) sg_qn_del(qt[t]); free(qt); }

/* Tier 4 over Q: Ntot+1 solves of perturbed models through ops->solve,
 * then the Lagrange value at t=0 in exact arithmetic.  If amp_bits is
 * given it receives the largest log2(sum_i |w_i g_i| / |g(0)|) over the
 * exported values: the number of bits by which an inexact solver's error
 * is amplified by the extrapolation. */
static inline int sg_interp_exact_amp(qnmodel* q, mpq_t* out, const sg_ops* ops, double* amp_bits)
{
	int nres = SG_NRES(q), D = sg_interp_degree(q), T = sg_interp_points(q);
	int got = 0, t, i, j, u, rc = 0;
	qnmodel** qt = sg_perturb_build(q, T);
	mpq_t** g = (mpq_t**) malloc((D+1)*sizeof(mpq_t*));
	int* ts = (int*) malloc((D+1)*sizeof(int));
	for (i=0;i<=D;i++) g[i] = sg_mpq_arr(nres);
	for (t=1; t<=T && got<D+1; t++)
	{
		int r1 = ops->solve(qt[t-1], g[got], ops->ctx);
		if (r1 == 0) ts[got++] = t;
		else if (r1 == 2) { rc = 2; break; }
	}
	if (rc == 0 && got < D+1) rc = 1;
	if (rc == 0)
	{
		mpq_t w, tmp, *mag = sg_mpq_arr(nres);
		mpq_init(w); mpq_init(tmp);
		for (u=0;u<nres;u++) mpq_set_ui(out[u], 0, 1);
		for (i=0;i<=D;i++)
		{
			mpq_set_ui(w, 1, 1);
			for (j=0;j<=D;j++) if (j != i)
			{
				long dd = (long)ts[i] - ts[j];           /* (0 - t_j)/(t_i - t_j) */
				mpq_set_si(tmp, dd < 0 ? ts[j] : -ts[j], (unsigned long)(dd < 0 ? -dd : dd));
				mpq_canonicalize(tmp);
				mpq_mul(w, w, tmp);
			}
			for (u=0;u<nres;u++)
			{
				mpq_mul(tmp, w, g[i][u]);
				mpq_add(out[u], out[u], tmp);
				if (amp_bits) { mpq_abs(tmp, tmp); mpq_add(mag[u], mag[u], tmp); }
			}
		}
		/* G(N-1_r) and G^{+k}(N-1_r) vanish identically when N_r = 0; an
		 * inexact solver would otherwise leave rounding noise there. */
		for (u=0;u<q->R;u++) if (q->N[u] == 0)
		{
			int kk;
			mpq_set_ui(out[1+u], 0, 1); mpq_set_ui(mag[1+u], 0, 1);
			for (kk=0;kk<q->M;kk++) { mpq_set_ui(out[SG_IDX_GK(q->R,kk,u)], 0, 1); mpq_set_ui(mag[SG_IDX_GK(q->R,kk,u)], 0, 1); }
		}
		if (amp_bits)
		{
			double a = 0.0;
			for (u=0;u<nres;u++)
			{
				double lm, lo;
				if (mpq_sgn(mag[u]) == 0) continue;
				lm = (double)mpz_sizeinbase(mpq_numref(mag[u]),2) - (double)mpz_sizeinbase(mpq_denref(mag[u]),2);
				if (mpq_sgn(out[u]) == 0) { a = 1e9; break; }   /* total cancellation */
				lo = (double)mpz_sizeinbase(mpq_numref(out[u]),2) - (double)mpz_sizeinbase(mpq_denref(out[u]),2);
				if (lm - lo + 2 > a) a = lm - lo + 2;
			}
			*amp_bits = a;
		}
		mpq_clear(w); mpq_clear(tmp); sg_mpq_arr_free(mag, nres);
	}
	for (i=0;i<=D;i++) sg_mpq_arr_free(g[i], nres);
	free(g); free(ts);
	sg_perturb_free(qt, T);
	return rc;
}

static inline int sg_interp_exact(qnmodel* q, mpq_t* out, void* vops)
{
	return sg_interp_exact_amp(q, out, (const sg_ops*) vops, NULL);
}

/* All tiers.  Returns 0 and the exact exported constants of q in out[],
 * 1 if every tier failed, or 2 on a hard error already reported.
 * skip_direct: the caller has already found the model as given singular. */
static inline int sg_solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz, const sg_ops* ops, int skip_direct)
{
	int rc;
	snprintf(tier, tsz, "direct");
	rc = sg_solve_reduce(q, out, tier, tsz, ops, skip_direct);
	if (rc != 1) return rc;
	rc = ops->interp ? ops->interp(q, out, ops->ctx) : sg_interp_exact(q, out, (void*) ops);
	if (rc == 0)
		snprintf(tier, tsz, "interpolation (%d+1 perturbed solves%s)", sg_interp_degree(q),
		         ops->interp_note ? ops->interp_note : "");
	return rc;
}

#endif
