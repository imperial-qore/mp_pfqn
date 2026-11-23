/*
 * replicas.h - the station multiplicity mi for the load-dependent solvers
 * (header-only).
 *
 * In a .qn file, mi is the number of identical copies of a station: mi
 * single-server queues, each with the station's demand L for every class
 * (and, with a MU section, each with the station's rates).  That is how
 * ca, mva, recal, the MoM/CoMoM family and the marginal solvers read it,
 * and the queue lengths they report for the station are the totals over
 * its copies.  A queue with c servers sharing one waiting line is written
 * with a MU section instead, rates min(n, c).
 *
 * gld, comomld and mvaldmx used to turn mi into the rates
 * min(n, mi) of one multiserver queue, a different model.  They now go
 * through one of the two routines below and solve the copies model.
 *
 * qn_replica_mu: mi copies of a station whose single copy has rates
 * mu(1..), a(t) = 1/(mu(1)...mu(t)), behave exactly as one load-dependent
 * station with
 *
 *   c = a * a * ... * a  (mi-fold convolution),  mu_agg(t) = c(t-1)/c(t).
 *
 * In a product-form state the copies hold class vectors n^1..n^mi summing
 * to n; their weights are prod_i (|n^i|!/prod_r n^i_r!) prod_r L_r^n^i_r
 * a(|n^i|), and summed over the splits with fixed totals t_i the
 * multinomials give |n|!/prod_r n_r! (arrange the |n| jobs in a row, cut
 * it into blocks of t_1, ..., t_mi).  So the station's weight is
 * (|n|!/prod_r n_r!) prod_r L_r^n_r c(|n|): a load-dependent station with
 * 1/(mu_agg(1)...mu_agg(t)) = c(t).  For load-independent copies
 * c(t) = C(t+mi-1, mi-1) and mu_agg(t) = t/(t+mi-1).  Exact; G, the
 * throughputs and the station totals are unchanged.
 *
 * The aggregated rates are not constant beyond Nt, so a solver that
 * extends them past Nt (mvaldmx, for open classes) uses qn_expand_replicas
 * instead: one station per copy, whose queue lengths are summed back.
 */
#ifndef REPLICAS_H
#define REPLICAS_H
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

/* Rates of every station under the copies meaning of mi, written to
 * qn->mu (allocated with unit rates when there is no MU section) and
 * qn->isLD set; mi is left as read.  Returns 0, or 1 if a station with
 * mi > 1 has a zero rate in its MU row. */
static inline int qn_replica_mu(qnmodel* qn)
{
	int M = qn->M, Nt = qn->Nt, m, t, s, i;
	if (Nt <= 0) return 0;
	if (qn->mu == NULL) {
		qn->mu = (mpq_t**) malloc((size_t)M * sizeof(mpq_t*));
		for (m = 0; m < M; m++) {
			qn->mu[m] = (mpq_t*) malloc((size_t)Nt * sizeof(mpq_t));
			for (t = 0; t < Nt; t++) mpq_init(qn->mu[m][t]), mpq_set_ui(qn->mu[m][t], 1, 1);
		}
	}
	qn->isLD = 1;
	mpq_t* a = (mpq_t*) malloc((size_t)(Nt + 1) * sizeof(mpq_t));
	mpq_t* c = (mpq_t*) malloc((size_t)(Nt + 1) * sizeof(mpq_t));
	mpq_t* w = (mpq_t*) malloc((size_t)(Nt + 1) * sizeof(mpq_t));
	mpq_t p;
	mpq_init(p);
	for (t = 0; t <= Nt; t++) { mpq_init(a[t]); mpq_init(c[t]); mpq_init(w[t]); }
	int rc = 0;
	for (m = 0; m < M && rc == 0; m++) {
		if (qn->mi[m] <= 1) continue;
		mpq_set_ui(a[0], 1, 1);
		for (t = 1; t <= Nt; t++) {
			if (mpq_sgn(qn->mu[m][t-1]) == 0) {
				fprintf(stderr, "station %d: a zero rate with multiplicity %d is not supported\n", m + 1, qn->mi[m]);
				rc = 1;
				break;
			}
			mpq_div(a[t], a[t-1], qn->mu[m][t-1]);
		}
		if (rc) break;
		/* c = a^{*mi} by repeated squaring: about 2 log2(mi) truncated
		 * convolutions instead of mi-1 (a ZRP with L=1024 sites is one
		 * station with mi=1024).  a is the running power a^{*2^j}, c
		 * collects the powers of the set bits of mi; both are exact, so
		 * the result is the same rational as the sequential product. */
		mpq_set_ui(c[0], 1, 1);
		for (t = 1; t <= Nt; t++) mpq_set_ui(c[t], 0, 1);
		for (i = qn->mi[m]; i > 0; i >>= 1) {
			if (i & 1) {
				for (t = 0; t <= Nt; t++) {
					mpq_set_ui(w[t], 0, 1);
					for (s = 0; s <= t; s++) { mpq_mul(p, c[s], a[t-s]); mpq_add(w[t], w[t], p); }
				}
				for (t = 0; t <= Nt; t++) mpq_swap(c[t], w[t]);
			}
			if (i > 1) {
				for (t = 0; t <= Nt; t++) {
					mpq_set_ui(w[t], 0, 1);
					for (s = 0; s <= t; s++) { mpq_mul(p, a[s], a[t-s]); mpq_add(w[t], w[t], p); }
				}
				for (t = 0; t <= Nt; t++) mpq_swap(a[t], w[t]);
			}
		}
		for (t = 1; t <= Nt; t++) mpq_div(qn->mu[m][t-1], c[t-1], c[t]);
	}
	for (t = 0; t <= Nt; t++) { mpq_clear(a[t]); mpq_clear(c[t]); mpq_clear(w[t]); }
	free(a); free(c); free(w);
	mpq_clear(p);
	return rc;
}

/* The copies model with one station per copy, mi = 1 each: demands and
 * (if any) rate rows are duplicated; N, Z and lambda are shared with qn.
 * *owner (length = new M) receives the original station of each copy.
 * Free with qn_free_expanded. */
static inline qnmodel* qn_expand_replicas(const qnmodel* qn, int** owner)
{
	int M = 0, m, i, r, t, k = 0;
	for (m = 0; m < qn->M; m++) M += qn->mi[m];
	qnmodel* e = (qnmodel*) calloc(1, sizeof(qnmodel));
	*e = *qn;
	e->M = M;
	e->L = (mpz_t**) malloc((size_t)M * sizeof(mpz_t*));
	e->mi = (int*) malloc((size_t)M * sizeof(int));
	e->mu = qn->mu ? (mpq_t**) malloc((size_t)M * sizeof(mpq_t*)) : NULL;
	*owner = (int*) malloc((size_t)M * sizeof(int));
	for (m = 0; m < qn->M; m++)
		for (i = 0; i < qn->mi[m]; i++, k++) {
			(*owner)[k] = m;
			e->mi[k] = 1;
			e->L[k] = (mpz_t*) malloc((size_t)qn->R * sizeof(mpz_t));
			for (r = 0; r < qn->R; r++) mpz_init_set(e->L[k][r], qn->L[m][r]);
			if (e->mu) {
				e->mu[k] = (mpq_t*) malloc((size_t)qn->Nt * sizeof(mpq_t));
				for (t = 0; t < qn->Nt; t++) { mpq_init(e->mu[k][t]); mpq_set(e->mu[k][t], qn->mu[m][t]); }
			}
		}
	return e;
}

/* Frees what qn_expand_replicas allocated (not the shared N, Z, lambda),
 * including a rate table a solver allocated on it. */
static inline void qn_free_expanded(qnmodel* e, int R)
{
	int k, r, t;
	for (k = 0; k < e->M; k++) {
		for (r = 0; r < R; r++) mpz_clear(e->L[k][r]);
		free(e->L[k]);
		if (e->mu) { for (t = 0; t < e->Nt; t++) mpq_clear(e->mu[k][t]); free(e->mu[k]); }
	}
	free(e->L); free(e->mi); free(e->mu); free(e);
}

#endif
