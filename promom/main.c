#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <gmpla.h>
#include "promom.h"

qnmodel* qnm;
double t0, t1;
struct rusage ruse;

/* ------------------------------------------------------------------ */

static mpq_vec_t* vecarray(int count, int len)
{
	int i;
	mpq_vec_t* v = (mpq_vec_t*) calloc(count, sizeof(mpq_vec_t));
	for (i = 0; i < count; i++) v[i] = mpq_vec(len, 0, 1);
	return v;
}

static void vecarray_free(mpq_vec_t* v, int count, int len)
{
	int i, j;
	if (!v) return;
	for (i = 0; i < count; i++) {
		for (j = 0; j < len + 1; j++) mpq_clear(v[i][j]);
		free(v[i]);
	}
	free(v);
}

/* dense rectangular mat-vec over a column range: rop[nrows] += mat * v */
static void matvec_range(mpq_vec_t rop, mpq_mat_t mat, int nrows, int c0, int c1, mpq_vec_t v)
{
	int i, j; mpq_t t; mpq_init(t);
	for (i = 0; i < nrows; i++)
		for (j = c0; j < c1; j++)
			if (mpq_sgn(mat[i][j]) && mpq_sgn(v[j])) {
				mpq_mul(t, mat[i][j], v[j]);
				mpq_add(rop[i], rop[i], t);
			}
	mpq_clear(t);
}

static void swap_stations(qnmodel* qn, int a, int b)
{
	if (a == b) return;
	mpz_t* tl = qn->L[a]; qn->L[a] = qn->L[b]; qn->L[b] = tl;
	int tm = qn->mi[a];   qn->mi[a] = qn->mi[b]; qn->mi[b] = tm;
}

/* copy a square submatrix mat[0..nr)[c0..c1) into a fresh dense matrix */
static mpq_mat_t submat(mpq_mat_t mat, int nr, int c0, int c1)
{
	int i, j;
	mpq_mat_t s = mpq_matzeros(nr, c1 - c0);
	for (i = 0; i < nr; i++)
		for (j = c0; j < c1; j++)
			mpq_set(s[i][j-c0], mat[i][j]);
	return s;
}

static void free_mat(mpq_mat_t m, int rows, int cols)
{
	int i, j;
	for (i = 0; i <= rows; i++) {
		for (j = 0; j <= cols; j++) mpq_clear(m[i][j]);
		free(m[i]);
	}
	free(m);
}

/* ------------------------------------------------------------------
 * solve_promom - the population sweep.
 *
 * Mirrors mom/main.c solve_model, but every basis vector carries the extra
 * index n = 0..sumN (jobs at the reference station), and every step system
 * gains the n-coupled blocks of genqmatrix.  Returns 0, or -1 on a singular
 * step matrix.
 */
static int solve_promom(qnmodel* qn, int sumN, mpq_t* dist)
{
	int M = qn->M, R = qn->R;
	int r, nr, n, i, s, t;
	int* Ncur = (int*) int_vec(R, 0);
	int ret = 0;

	mpq_vec_t* g   = NULL;   /* basis at the current point, per n     */
	mpq_vec_t* gr  = NULL;   /* basis at the previous point, per n    */
	int len = 0;

	for (r = 1; r <= R; r++) {
		qcombs* Ik = qcombs_new(M, r);
		qcombs* I  = qcombs_new(M, r-1);
		int cardk = Ik->card, cardi = I->card;
		int off = cardk * r;
		int newlen = (cardk + cardi) * r;

		if (r == 1) {
			/* population 0: every constant of the empty model is 1, and
			 * the reference station holds no jobs, so only n = 0 is set */
			g  = vecarray(sumN+1, newlen);
			gr = vecarray(sumN+1, newlen);
			for (i = 0; i < newlen; i++) mpq_set_si(g[0][i], 1, 1);
			len = newlen;
		} else {
			/* class transition: the old top level becomes the new
			 * level-(r-1) block; its class-(r-1) decrement component comes
			 * from the basis one population step back (gr).  The new
			 * level-r block then follows from the CE and class-s PC rows,
			 * which at class-r population 0 have no previous-point term. */
			QMatrices* qm = genqmatrix(Ik, I, qn, Ncur, r);
			int nr12 = off;                     /* CE + class-s PC rows */
			mpq_mat_t A11 = submat(qm->A, nr12, 0, off);
			int* lu = mpq_ludcmp(A11, off);
			if (!lu) {
				fprintf(stderr, "promom: singular transition matrix at class %d\n", r);
				free_mat(A11, off, off); free_qmatrices(qm);
				qcombs_free(Ik); qcombs_free(I); ret = -1; goto done;
			}

			mpq_vec_t* gnew = vecarray(sumN+1, newlen);
			for (n = 0; n <= sumN; n++) {
				/* known level-(r-1) block */
				for (i = 0; i < cardi; i++) {
					for (s = 0; s <= r-2; s++)
						mpq_set(gnew[n][off + i*r + s], g[n][i*(r-1) + s]);
					mpq_set(gnew[n][off + i*r + (r-1)], gr[n][i*(r-1)]);
				}
				/* rhs = n * DC * gnew(n-1) - A12 * (level-(r-1) block) */
				mpq_vec_t rhs = mpq_vec(off, 0, 1);
				if (n > 0) {
					mpq_vec_t acc = mpq_vec(off, 0, 1);
					matvec_range(acc, qm->DC, nr12, off, newlen, gnew[n-1]);
					mpq_t nq; mpq_init(nq); mpq_set_si(nq, n, 1);
					for (t = 0; t < off; t++) { mpq_mul(acc[t], acc[t], nq); mpq_add(rhs[t], rhs[t], acc[t]); }
					mpq_clear(nq);
					for (t = 0; t < off+1; t++) mpq_clear(acc[t]);
					free(acc);
				}
				{
					mpq_vec_t acc = mpq_vec(off, 0, 1);
					matvec_range(acc, qm->A, nr12, off, newlen, gnew[n]);
					for (t = 0; t < off; t++) mpq_sub(rhs[t], rhs[t], acc[t]);
					for (t = 0; t < off+1; t++) mpq_clear(acc[t]);
					free(acc);
				}
				if (mpq_lubksb(A11, rhs, off, lu) < 0) {
					fprintf(stderr, "promom: singular transition back-substitution\n");
					ret = -1;
				}
				for (t = 0; t < off; t++) mpq_set(gnew[n][t], rhs[t]);
				for (t = 0; t < off+1; t++) mpq_clear(rhs[t]);
				free(rhs);
			}
			free(lu); free_mat(A11, off, off); free_qmatrices(qm);
			vecarray_free(g, sumN+1, len);
			vecarray_free(gr, sumN+1, len);
			g  = gnew;
			gr = vecarray(sumN+1, newlen);
			len = newlen;
			if (ret) { qcombs_free(Ik); qcombs_free(I); goto done; }
		}

		for (nr = 1; nr <= qn->N[r-1]; nr++) {
			Ncur[r-1] = nr;
			int sumNcur = sum(Ncur, R);

			/* the previous point is the current one before this step */
			for (n = 0; n <= sumN; n++) mpq_vecdup(gr[n], g[n], len);

			QMatrices* qm = genqmatrix(Ik, I, qn, Ncur, r);
			if (qm->nrows != qm->ncols) {
				fprintf(stderr, "promom: system not square (%d x %d) at class %d\n",
				        qm->nrows, qm->ncols, r);
				free_qmatrices(qm); qcombs_free(Ik); qcombs_free(I); ret = -1; goto done;
			}
			mpq_mat_t A = submat(qm->A, qm->nrows, 0, qm->ncols);
			int* lu = mpq_ludcmp(A, qm->ncols);
			if (!lu) {
				fprintf(stderr, "promom: singular step matrix at class %d, Nr=%d\n", r, nr);
				free_mat(A, qm->nrows, qm->ncols); free_qmatrices(qm);
				qcombs_free(Ik); qcombs_free(I); ret = -1; goto done;
			}

			for (n = 0; n <= sumN; n++) {
				mpq_vec_t rhs = mpq_vec(len, 0, 1);
				if (n <= sumNcur) {
					matvec_range(rhs, qm->B, qm->nrows, 0, len, gr[n]);
					if (n > 0) {
						mpq_vec_t acc = mpq_vec(len, 0, 1);
						matvec_range(acc, qm->DC, qm->nrows, 0, len, g[n-1]);
						matvec_range(acc, qm->DD, qm->nrows, 0, len, gr[n-1]);
						mpq_t nq; mpq_init(nq); mpq_set_si(nq, n, 1);
						for (t = 0; t < len; t++) { mpq_mul(acc[t], acc[t], nq); mpq_add(rhs[t], rhs[t], acc[t]); }
						mpq_clear(nq);
						for (t = 0; t < len+1; t++) mpq_clear(acc[t]);
						free(acc);
					}
					if (mpq_lubksb(A, rhs, len, lu) < 0) {
						fprintf(stderr, "promom: singular back-substitution\n");
						ret = -1;
					}
				}
				for (t = 0; t < len; t++) mpq_set(g[n][t], rhs[t]);
				for (t = 0; t < len+1; t++) mpq_clear(rhs[t]);
				free(rhs);
				if (ret) break;
			}
			free(lu); free_mat(A, qm->nrows, qm->ncols); free_qmatrices(qm);
			if (ret) { qcombs_free(Ik); qcombs_free(I); goto done; }
		}
		qcombs_free(Ik); qcombs_free(I);
	}

	ret = qdecrease(qn, g, g, gr, sumN, dist);

done:
	vecarray_free(g,  sumN+1, len);
	vecarray_free(gr, sumN+1, len);
	free(Ncur);
	return ret;
}

/* ------------------------------------------------------------------ */

int main(int argc, char** argv)
{
	t0 = CPUTIME;
	int i;
	bool queue_output = false, prob_output = false;
	char* model_file = NULL;

	for (i = 1; i < argc; i++) {
		if      (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) queue_output = true;
		else if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--prob")) prob_output = true;
		else if (argv[i][0] != '-') model_file = argv[i];
	}
	if (!model_file) {
		printf("USAGE: %s [-q|--qlen] [-P|--prob] model.qn\n", argv[0]);
		printf("  promom: marginal queue-length probabilities on the MoM basis,\n");
		printf("  the MoM-basis counterpart of bin/procomom.\n");
		printf("  -q, --qlen : mean queue lengths, one per station\n");
		printf("  -P, --prob : marginal distribution per station, one row per station\n");
		return -1;
	}

	qnm = (qnmodel*) readmodel(model_file);
	int M = qnm->M, R = qnm->R;
	{ int c; for (c = 0; c < R; c++) if (qnm->N[c] < 0) {
		fprintf(stderr, "promom: open/mixed models are not supported\n"); return 1; } }
	if (qnm->hasOpen) { fprintf(stderr, "promom: open/mixed models are not supported\n"); return 1; }
	if (qnm->isLD)    { fprintf(stderr, "promom: load-dependent stations are not supported\n"); return 1; }

	nckinit(M + R + 1, M);
	int sumN = sum(qnm->N, R);

	mpq_t* dist = (mpq_t*) calloc(sumN+1, sizeof(mpq_t));
	for (i = 0; i <= sumN; i++) mpq_init(dist[i]);

	mpf_t fval; mpf_init(fval);
	int k, j;
	if (!queue_output && !prob_output)
		printf("========== Marginal Queue-Length Probabilities ==========\n");

	for (k = 1; k <= M; k++) {
		swap_stations(qnm, k-1, M-1);
		int rc = solve_promom(qnm, sumN, dist);
		swap_stations(qnm, k-1, M-1);
		if (rc < 0) {
			fprintf(stderr, "promom: singular system; no exact answer for station %d\n", k);
			return 1;
		}

		mpq_t total, Qk, term, p;
		mpq_init(total); mpq_init(Qk); mpq_init(term); mpq_init(p);
		mpq_set_si(total, 0, 1); mpq_set_si(Qk, 0, 1);
		for (j = 0; j <= sumN; j++) mpq_add(total, total, dist[j]);

		if (!queue_output && !prob_output) printf("\nStation %d:\n", k);
		for (j = 0; j <= sumN; j++) {
			mpq_set(p, dist[j]);
			if (mpq_sgn(total) != 0) mpq_div(p, p, total);
			if (mpq_sgn(p) != 0) {
				mpq_set_si(term, j, 1);
				mpq_mul(term, term, p);
				mpq_add(Qk, Qk, term);
			}
			mpf_set_q(fval, p);
			if (prob_output) {
				printf("%.15e", mpf_get_d(fval));
				if (j < sumN) printf(" ");
			} else if (!queue_output) {
				double pv = mpf_get_d(fval);
				if (pv > 1e-20 || j <= 5) printf("  P(n_%d = %d) = %.15e\n", k, j, pv);
			}
		}
		if (prob_output) printf("\n");
		/* mi_k > 1 means mi_k replicated single-server queues (this is how
		 * bin/ca expands them), and the distribution above is the marginal
		 * at ONE copy.  The station total is therefore mi_k times its mean,
		 * which is what -q reports so it lines up with ca/comom. */
		{ mpq_t mq; mpq_init(mq); mpq_set_si(mq, qnm->mi[k-1], 1);
		  mpq_mul(Qk, Qk, mq); mpq_clear(mq); }
		mpf_set_q(fval, Qk);
		if (queue_output)        printf("%.15e\n", mpf_get_d(fval));
		else if (!prob_output)   printf("  Q[%d] = %.15e\n", k, mpf_get_d(fval));

		mpq_clear(total); mpq_clear(Qk); mpq_clear(term); mpq_clear(p);
	}
	mpf_clear(fval);

	if (!queue_output && !prob_output) {
		printf("\n=========================================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (ProMoM): %g s\n", t1 - t0);
	}
	for (i = 0; i <= sumN; i++) mpq_clear(dist[i]);
	free(dist);
	return 0;
}
