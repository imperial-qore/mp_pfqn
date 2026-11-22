#ifndef SINGULAR_MARG_H
#define SINGULAR_MARG_H

/* Exact fallback for the marginal-probability solvers (promom, procomom)
 * on models where their linear systems are singular.
 *
 * The counterpart of singular.h, whose model transformations it reuses,
 * for an output that is a distribution per station rather than a set of
 * normalizing constants.  The solver supplies a callback that fills
 *   out[k*(Nt+1) + n] = q_k(N, n),   k = 0..M-1, n = 0..Nt,
 * the UNNORMALIZED marginal at one copy of station k (one row per
 * station, Nt = total population), or reports that its system is
 * singular.  sgm_solve_tiered() returns the same table NORMALIZED (each
 * row sums to 1), trying cheapest first:
 *   1. other class orders: the marginal at a station carries no class
 *      index, so the result is used as it comes back;
 *   2. stations with identical demand rows folded into one station of
 *      summed multiplicity: every copy of the folded station is a copy of
 *      each member, so member k reads row g(k);
 *   3. classes with proportional demands and think times aggregated: the
 *      generating function of q_k(N, n) depends on such a group only
 *      through x_s + sum_i c_i x_i (as for G, see sg_merge_classes), so
 *      q_k(N, n) = q'_k(n_agg, n) F(N) and the normalized row is the same;
 *   4. the exact limit of a perturbation.  q_k(N, n) is a homogeneous
 *      polynomial of degree Nt in (L, Z) jointly, like G, so along
 *      (L, Z) + t (DL, DZ) it is a polynomial of degree <= Nt in t, and
 *      the Lagrange interpolant through Nt+1 successful perturbed solves
 *      evaluated at t = 0 gives it exactly.  The normalized probabilities
 *      are rational in t, not polynomial, which is why the callback must
 *      return the unnormalized values.  One further solve is predicted
 *      from the interpolant and compared exactly: if a solver's scaling
 *      ever made its output non-polynomial, the tier reports failure
 *      instead of returning a wrong limit.
 * Tiers 1-3 are applied again to each reduced model.  Every tier is
 * exact. */

#include "singular.h"

#define SGM_NRES(q) ((q)->M * ((q)->Nt + 1))

typedef struct {
	sg_solve_fn solve;   /* same contract as singular.h: 0 / 1 singular / 2 error */
	void* ctx;
} sgm_ops;

static inline void sgm_fix_nt(qnmodel* q) { sg_qn_fix_nt(q); }

/* Tier 1: the model as given (unless skip_direct), then other class orders */
static inline int sgm_solve_order(qnmodel* q, mpq_t* out, int* reordered, const sgm_ops* ops, int skip_direct)
{
	int R = q->R, M = q->M, i, j, k, tries = 1, rc = 1;
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
		for (j=0;j<R;j++) {
			p->N[j] = q->N[perm[j]];
			mpz_set(p->Z[j], q->Z[perm[j]]);
			for (k=0;k<M;k++) mpz_set(p->L[k][j], q->L[k][perm[j]]);
		}
		for (k=0;k<M;k++) p->mi[k] = q->mi[k];
		sg_qn_fix_nt(p);
		rc = ops->solve(p, out, ops->ctx);
		if (rc == 0) *reordered = 1;
		sg_qn_del(p);
	}
	free(perm);
	return rc;
}

/* Tiers 1-3; out is unnormalized up to one scale factor per station */
static inline int sgm_solve_reduce(qnmodel* q, mpq_t* out, char* tier, size_t tsz, const sgm_ops* ops, int skip_direct)
{
	int R = q->R, M = q->M, W = q->Nt + 1, rc, reordered, k, n;
	rc = sgm_solve_order(q, out, &reordered, ops, skip_direct);
	if (rc == 0) { if (reordered) snprintf(tier, tsz, "class order"); return 0; }
	if (rc != 1) return rc;

	/* identical stations */
	{
		int* grp = (int*) malloc(M*sizeof(int));
		qnmodel* ms = sg_merge_stations(q, grp);
		if (ms)
		{
			int nres2 = SGM_NRES(ms);
			mpq_t* o = sg_mpq_arr(nres2);
			char sub[96] = "";
			rc = sgm_solve_reduce(ms, o, sub, sizeof sub, ops, 0);
			if (rc == 0)
			{
				for (k=0;k<M;k++) for (n=0;n<W;n++) mpq_set(out[k*W + n], o[grp[k]*W + n]);
				snprintf(tier, tsz, "merged stations (M=%d->%d)%s%s", M, ms->M, sub[0] ? ", " : "", sub);
			}
			sg_mpq_arr_free(o, nres2);
			sg_qn_del(ms);
			if (rc != 1) { free(grp); return rc; }
		}
		free(grp);
	}

	/* proportional classes: same Nt and M, rows differ by F(N) only */
	{
		int* grp = (int*) malloc(R*sizeof(int));
		mpq_t* c = sg_mpq_arr(R);
		qnmodel* mc = sg_merge_classes(q, grp, c);
		if (mc)
		{
			char sub[96] = "";
			rc = sgm_solve_reduce(mc, out, sub, sizeof sub, ops, 0);
			if (rc == 0)
				snprintf(tier, tsz, "merged classes (R=%d->%d)%s%s", R, mc->R, sub[0] ? ", " : "", sub);
			sg_qn_del(mc);
		}
		sg_mpq_arr_free(c, R); free(grp);
	}
	return rc;
}

/* value at x of the Lagrange interpolant through (ts[i], g[i][u]), i<np */
static inline void sgm_lagrange(mpq_t* res, int nres, mpq_t** g, const int* ts, int np, long x)
{
	int i, j, u;
	mpq_t w, tmp; mpq_init(w); mpq_init(tmp);
	for (u=0;u<nres;u++) mpq_set_ui(res[u], 0, 1);
	for (i=0;i<np;i++)
	{
		mpq_set_ui(w, 1, 1);
		for (j=0;j<np;j++) if (j != i)
		{
			long num = x - ts[j], den = (long)ts[i] - ts[j];
			if (den < 0) { num = -num; den = -den; }
			mpq_set_si(tmp, num, (unsigned long)den);
			mpq_canonicalize(tmp);
			mpq_mul(w, w, tmp);
		}
		for (u=0;u<nres;u++) if (mpq_sgn(g[i][u])) { mpq_mul(tmp, w, g[i][u]); mpq_add(res[u], res[u], tmp); }
	}
	mpq_clear(w); mpq_clear(tmp);
}

/* Tier 4: Nt+1 perturbed solves, one more as a check, then the value at 0 */
static inline int sgm_interp(qnmodel* q, mpq_t* out, const sgm_ops* ops)
{
	int nres = SGM_NRES(q), D = sg_interp_degree(q), T = sg_interp_points(q) + 1;
	int need = D + 2, got = 0, t, i, u, rc = 0;
	qnmodel** qt = sg_perturb_build(q, T);
	mpq_t** g = (mpq_t**) malloc(need*sizeof(mpq_t*));
	int* ts = (int*) malloc(need*sizeof(int));
	for (i=0;i<need;i++) g[i] = sg_mpq_arr(nres);
	for (t=1; t<=T && got<need; t++)
	{
		int r1 = ops->solve(qt[t-1], g[got], ops->ctx);
		if (r1 == 0) ts[got++] = t;
		else if (r1 == 2) { rc = 2; break; }
	}
	if (rc == 0 && got < need) rc = 1;
	if (rc == 0)
	{
		mpq_t* chk = sg_mpq_arr(nres);
		sgm_lagrange(chk, nres, g, ts, D+1, ts[D+1]);
		for (u=0;u<nres;u++) if (!mpq_equal(chk[u], g[D+1][u])) break;
		if (u < nres) {
			fprintf(stderr, "singular_marg: perturbed marginals are not polynomial of degree %d; "
			                "interpolation tier skipped\n", D);
			rc = 1;
		} else
			sgm_lagrange(out, nres, g, ts, D+1, 0);
		sg_mpq_arr_free(chk, nres);
	}
	for (i=0;i<need;i++) sg_mpq_arr_free(g[i], nres);
	free(g); free(ts);
	sg_perturb_free(qt, T);
	return rc;
}

/* each row of out divided by its sum (rows that sum to 0 are left as is) */
static inline void sgm_normalize(const qnmodel* q, mpq_t* out)
{
	int W = q->Nt + 1, k, n;
	mpq_t s; mpq_init(s);
	for (k=0;k<q->M;k++)
	{
		mpq_set_ui(s, 0, 1);
		for (n=0;n<W;n++) mpq_add(s, s, out[k*W + n]);
		if (mpq_sgn(s)) for (n=0;n<W;n++) mpq_div(out[k*W + n], out[k*W + n], s);
	}
	mpq_clear(s);
}

/* All tiers.  Returns 0 and the normalized per-copy marginals of q in
 * out[SGM_NRES(q)], 1 if every tier failed, 2 on a hard error already
 * reported.  skip_direct: the caller has already found q singular. */
static inline int sgm_solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz, const sgm_ops* ops, int skip_direct)
{
	int rc;
	sgm_fix_nt(q);
	snprintf(tier, tsz, "direct");
	rc = sgm_solve_reduce(q, out, tier, tsz, ops, skip_direct);
	if (rc == 1) {
		rc = sgm_interp(q, out, ops);
		if (rc == 0)
			snprintf(tier, tsz, "interpolation (%d+1 perturbed solves, 1 check)", sg_interp_degree(q));
	}
	if (rc == 0) sgm_normalize(q, out);
	return rc;
}

#endif
