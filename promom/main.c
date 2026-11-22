#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <gmpla.h>
#include "promom.h"
#include "singular_marg.h"

qnmodel* qnm;
double t0, t1;
struct rusage ruse;
/* set once the direct solve has failed: the fallback tiers make many
 * solves that may each be singular, and only their outcome is reported */
static int promom_quiet = 0;

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
 * Lower replica levels for mixed think times.
 *
 * qdecrease recovers a class-s decrement with Z_s = 0 from the conditioned
 * server-count identity, which is the sum of the PC over all classes and
 * so holds only when every Z is zero.  With mixed think times the summed
 * PC keeps the terms Z_u q(m, Q-1_u, n) at two decrements below N, which
 * the top two levels do not carry.  As mom does since 14a8a4d, the sweep
 * then carries levels 0..r-2 as well, each advanced along class r by the
 * class-r PC (the rows-3 identity of genqmatrix at a lower level),
 *
 *   nr q_l(i, P-1_s, n) = Z_r q_l(i, P-1_s-1_r, n)
 *       + sum_j m'_j L_jr q_{l+1}(i+1_j, P-1_s-1_r, n)
 *       + n L_Mr q_l(i, P-1_s-1_r, n-1),
 *
 * every right-hand term at the previous point, and level 0 at N is read
 * directly, so no descent is needed.  low[l][n] holds level l with r
 * components per combination, laid out as the level-(r-1) block of g.  */
typedef struct {
	int nl;          /* levels 0..nl-1 (nl = R-1)                         */
	int r;           /* components per combination of the stored vectors */
	qcombs** C;      /* combinations at levels 0..nl (nl = top lookup)    */
	int** up;        /* up[l][i*M+j]: index of C[l] row i plus 1_j in C[l+1] */
	mpq_vec_t** low; /* low[l][n], current point                         */
	mpq_vec_t** lowr;/* low[l][n], previous point                        */
} lowlv;

static lowlv* lowlv_new(int M, int R)
{
	int l, i, j, t;
	lowlv* lv = (lowlv*) calloc(1, sizeof(lowlv));
	lv->nl = R - 1;
	lv->C  = (qcombs**) calloc(R, sizeof(qcombs*));
	lv->up = (int**) calloc(R, sizeof(int*));
	lv->low  = (mpq_vec_t**) calloc(R, sizeof(mpq_vec_t*));
	lv->lowr = (mpq_vec_t**) calloc(R, sizeof(mpq_vec_t*));
	int* sh = (int*) calloc(M, sizeof(int));
	for (l = 0; l <= lv->nl; l++) lv->C[l] = qcombs_new(M, l);
	for (l = 0; l < lv->nl; l++) {
		lv->up[l] = (int*) calloc(lv->C[l]->card * M, sizeof(int));
		for (i = 0; i < lv->C[l]->card; i++)
			for (j = 0; j < M; j++) {
				for (t = 0; t < M; t++) sh[t] = lv->C[l]->combs[i][t];
				sh[j]++;
				lv->up[l][i*M + j] = qcombs_index(lv->C[l+1], sh);
			}
	}
	free(sh);
	return lv;
}

static void lowlv_free(lowlv* lv, int sumN)
{
	int l;
	if (!lv) return;
	for (l = 0; l < lv->nl; l++) {
		if (lv->low[l])  vecarray_free(lv->low[l],  sumN+1, lv->C[l]->card * lv->r);
		if (lv->lowr[l]) vecarray_free(lv->lowr[l], sumN+1, lv->C[l]->card * lv->r);
		free(lv->up[l]);
	}
	for (l = 0; l <= lv->nl; l++) qcombs_free(lv->C[l]);
	free(lv->C); free(lv->up); free(lv->low); free(lv->lowr);
	free(lv);
}

/* class transition r-1 -> r: levels 0..r-2 gain the class-(r-1) decrement
 * from the previous-point snapshot.  Level r-2 comes from the second part
 * of the old g/gr (offset cardi*(r-1), cardi = card of level r-1). */
static void lowlv_seed(lowlv* lv, int r, int sumN, mpq_vec_t* g, mpq_vec_t* gr, int cardi)
{
	int l, n, i, s;
	for (l = 0; l <= r-2; l++) {
		int card = lv->C[l]->card;
		mpq_vec_t* nw = vecarray(sumN+1, card * r);
		for (n = 0; n <= sumN; n++) {
			mpq_t* old  = (l == r-2) ? g[n]  + cardi*(r-1) : lv->low[l][n];
			mpq_t* snap = (l == r-2) ? gr[n] + cardi*(r-1) : lv->lowr[l][n];
			for (i = 0; i < card; i++) {
				for (s = 0; s <= r-2; s++) mpq_set(nw[n][i*r + s], old[i*(r-1) + s]);
				mpq_set(nw[n][i*r + r-1], snap[i*(r-1)]);
			}
		}
		if (l < r-2) {
			vecarray_free(lv->low[l],  sumN+1, card * (r-1));
			vecarray_free(lv->lowr[l], sumN+1, card * (r-1));
		}
		lv->low[l]  = nw;
		lv->lowr[l] = vecarray(sumN+1, card * r);   /* stays 0 if N_r = 0 */
	}
	lv->r = r;
}

/* one population step of class r to nr: snapshot, then the class-r PC.
 * gr is the previous-point basis, its level r-1 block at offset off. */
static void lowlv_step(lowlv* lv, qnmodel* qn, int r, int nr, int sumN, mpq_vec_t* gr, int off)
{
	int M = qn->M, l, n, i, s, j;
	mpq_t acc, c, w;
	mpq_init(acc); mpq_init(c); mpq_init(w);
	for (l = 0; l <= r-2; l++)
		for (n = 0; n <= sumN; n++)
			mpq_vecdup(lv->lowr[l][n], lv->low[l][n], lv->C[l]->card * r);
	for (l = 0; l <= r-2; l++) {
		int card = lv->C[l]->card;
		for (n = 0; n <= sumN; n++) {
			mpq_t* above = (l+1 == r-1) ? gr[n] + off : lv->lowr[l+1][n];
			for (i = 0; i < card; i++) {
				for (s = 0; s < r; s++) {
					mpq_set_z(c, qn->Z[r-1]);
					mpq_mul(acc, c, lv->lowr[l][n][i*r + s]);
					for (j = 0; j < M; j++) {
						int mult = qn->mi[j] + lv->C[l]->combs[i][j] - (j == M-1 ? 1 : 0);
						if (mult == 0 || mpz_sgn(qn->L[j][r-1]) == 0) continue;
						mpq_set_z(c, qn->L[j][r-1]);
						mpq_set_si(w, mult, 1);
						mpq_mul(c, c, w);
						mpq_mul(w, c, above[lv->up[l][i*M + j]*r + s]);
						mpq_add(acc, acc, w);
					}
					if (n > 0 && mpz_sgn(qn->L[M-1][r-1]) != 0) {
						mpq_set_z(c, qn->L[M-1][r-1]);
						mpq_set_si(w, n, 1);
						mpq_mul(c, c, w);
						mpq_mul(w, c, lv->lowr[l][n-1][i*r + s]);
						mpq_add(acc, acc, w);
					}
					mpq_set_si(w, nr, 1);
					mpq_div(lv->low[l][n][i*r + s], acc, w);
				}
			}
		}
	}
	mpq_clear(acc); mpq_clear(c); mpq_clear(w);
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

	/* mixed think times: carry the lower levels (see lowlv above) */
	lowlv* lv = NULL;
	{
		int z0 = 0, zp = 0;
		for (s = 0; s < R; s++) { if (mpz_sgn(qn->Z[s]) == 0) z0 = 1; else zp = 1; }
		if (z0 && zp) lv = lowlv_new(M, R);
	}

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
				if (!promom_quiet) fprintf(stderr, "promom: singular transition matrix at class %d\n", r);
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
					if (!promom_quiet) fprintf(stderr, "promom: singular transition back-substitution\n");
					ret = -1;
				}
				for (t = 0; t < off; t++) mpq_set(gnew[n][t], rhs[t]);
				for (t = 0; t < off+1; t++) mpq_clear(rhs[t]);
				free(rhs);
			}
			free(lu); free_mat(A11, off, off); free_qmatrices(qm);
			if (lv) lowlv_seed(lv, r, sumN, g, gr, cardi);
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
			if (lv) lowlv_step(lv, qn, r, nr, sumN, gr, off);

			QMatrices* qm = genqmatrix(Ik, I, qn, Ncur, r);
			if (qm->nrows != qm->ncols) {
				if (!promom_quiet) fprintf(stderr, "promom: system not square (%d x %d) at class %d\n",
				        qm->nrows, qm->ncols, r);
				free_qmatrices(qm); qcombs_free(Ik); qcombs_free(I); ret = -1; goto done;
			}
			mpq_mat_t A = submat(qm->A, qm->nrows, 0, qm->ncols);
			int* lu = mpq_ludcmp(A, qm->ncols);
			if (!lu) {
				if (!promom_quiet) fprintf(stderr, "promom: singular step matrix at class %d, Nr=%d\n", r, nr);
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
						if (!promom_quiet) fprintf(stderr, "promom: singular back-substitution\n");
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

	if (lv) {
		/* level 0 at N, base component: q(N, n); lv exists only for R >= 2 */
		for (n = 0; n <= sumN; n++) mpq_set(dist[n], lv->low[0][n][0]);
		ret = 0;
	} else
		ret = qdecrease(qn, g, g, gr, sumN, dist);

done:
	lowlv_free(lv, sumN);
	vecarray_free(g,  sumN+1, len);
	vecarray_free(gr, sumN+1, len);
	free(Ncur);
	return ret;
}

/* ------------------------------------------------------------------
 * Every station's unnormalized marginal, out[k*(sumN+1)+n]: station k is
 * rotated into the reference position M.  The callback of singular_marg.h:
 * 0 solved, 1 singular. */
static int promom_all(qnmodel* qn, mpq_t* out, void* ctx)
{
	int M = qn->M, W = sum(qn->N, qn->R) + 1, k, n, rc = 0;
	(void) ctx;
	mpq_t* dist = sg_mpq_arr(W);
	for (k = 0; k < M && rc == 0; k++) {
		swap_stations(qn, k, M-1);
		rc = solve_promom(qn, W-1, dist) < 0 ? 1 : 0;
		swap_stations(qn, k, M-1);
		if (rc == 0) for (n = 0; n < W; n++) mpq_set(out[k*W + n], dist[n]);
	}
	sg_mpq_arr_free(dist, W);
	return rc;
}

int main(int argc, char** argv)
{
	t0 = CPUTIME;
	int i;
	bool queue_output = false, prob_output = false, no_fallback = false;
	char* model_file = NULL;

	for (i = 1; i < argc; i++) {
		if      (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) queue_output = true;
		else if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--prob")) prob_output = true;
		else if (!strcmp(argv[i], "-X") || !strcmp(argv[i], "--no-perturb")) no_fallback = true;
		else if (argv[i][0] != '-') model_file = argv[i];
	}
	if (!model_file) {
		printf("USAGE: %s [-q|--qlen] [-P|--prob] [-X] model.qn\n", argv[0]);
		printf("  promom: marginal queue-length probabilities on the MoM basis,\n");
		printf("  the MoM-basis counterpart of bin/procomom.\n");
		printf("  -q, --qlen : mean queue lengths, one per station\n");
		printf("  -P, --prob : marginal distribution per station, one row per station\n");
		printf("  -X, --no-perturb : fail on a singular system instead of running the\n");
		printf("               exact fallback (class orders, merged stations/classes,\n");
		printf("               interpolation; util/singular_marg.h)\n");
		return -1;
	}

	qnm = (qnmodel*) readmodel(model_file);
	int M = qnm->M, R = qnm->R;
	{ int c; for (c = 0; c < R; c++) if (qnm->N[c] < 0) {
		fprintf(stderr, "promom: open/mixed models are not supported\n"); return 1; } }
	if (qnm->hasOpen) { fprintf(stderr, "promom: open/mixed models are not supported\n"); return 1; }
	if (qnm->isLD)    { fprintf(stderr, "promom: load-dependent stations are not supported\n"); return 1; }

	nckinit(M + R + 1, M);
	int sumN = sum(qnm->N, R), W = sumN + 1;
	sg_qn_fix_nt(qnm);

	/* all stations first, so a singular system is known before printing */
	mpq_t* tab = sg_mpq_arr(M * W);
	char tier[256] = "direct";
	int rc = promom_all(qnm, tab, NULL);
	if (rc == 0)
		sgm_normalize(qnm, tab);
	else if (no_fallback) {
		fprintf(stderr, "promom: singular system; no exact answer (-X suppresses the exact fallback)\n");
		return 1;
	} else {
		sgm_ops ops = { promom_all, NULL };
		promom_quiet = 1;
		rc = sgm_solve_tiered(qnm, tab, tier, sizeof tier, &ops, 1);
		if (rc != 0) {
			fprintf(stderr, "promom: singular system, and every exact fallback failed "
			                "(class orders, merged stations or classes, interpolation)\n");
			return 1;
		}
		fprintf(stderr, "promom: singular system, solved exactly via %s\n", tier);
	}

	mpf_t fval; mpf_init(fval);
	int k, j;
	if (!queue_output && !prob_output) {
		printf("========== Marginal Queue-Length Probabilities ==========\n");
		if (strcmp(tier, "direct")) printf("\nSingular system, solved exactly via %s\n", tier);
	}

	for (k = 1; k <= M; k++) {
		mpq_t Qk, term;
		mpq_init(Qk); mpq_init(term);
		mpq_t* p = tab + (k-1)*W;

		if (!queue_output && !prob_output) printf("\nStation %d:\n", k);
		for (j = 0; j <= sumN; j++) {
			if (mpq_sgn(p[j]) != 0) {
				mpq_set_si(term, j, 1);
				mpq_mul(term, term, p[j]);
				mpq_add(Qk, Qk, term);
			}
			mpf_set_q(fval, p[j]);
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

		mpq_clear(Qk); mpq_clear(term);
	}
	mpf_clear(fval);

	if (!queue_output && !prob_output) {
		printf("\n=========================================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (ProMoM): %g s\n", t1 - t0);
	}
	sg_mpq_arr_free(tab, M * W);
	return 0;
}
