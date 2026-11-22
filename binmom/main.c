/*
 * binmom - MoM with the population steps of each class multiplied together
 * by binary splitting.
 *
 * Within class r the MoM coefficient blocks A11, A12, B1r and B2r do not
 * depend on n_r (bin/mom factors them once per class), and one step is
 *
 *   G(n)  = B2r g(n-1) / n,     Gk(n) = A11^{-1} (B1r g(n-1) - A12 G(n)),
 *
 * so with g = [Gk; G]
 *
 *   n g(n) = (n P + Q) g(n-1),  P = [A11^{-1} B1r; 0],
 *                               Q = [-A11^{-1} A12 B2r; B2r].
 *
 * Clearing the common denominator D of P and Q gives integer matrices
 * P' = D P, Q' = D Q and
 *
 *   g(N_r) = (N_r P' + Q') ... (P' + Q') g(0) / (D^N_r N_r!).
 *
 * bin/mom applies the N_r steps one at a time to a vector whose entries
 * grow by Theta(log n) bits per step, about d^2 N^2 log N bit operations
 * per class (d = basis size).  Here the product of the integer matrices is
 * formed by binary splitting, so the big multiplications happen between
 * operands of balanced size and GMP's subquadratic multiplication applies:
 * about d^omega M(N log N) bit operations.  It wins when N_r is large
 * against d and loses when d is large against N_r: the products are dense
 * d x d with entries widened by D, while mom's steps use sparse blocks of
 * about d*k entries, so the break-even population grows like d^2.
 * Measured on single-server models it lies at N_r ~ 1.9 d^2 (d = 28),
 * 2.2 d^2 (d = 48) and 2.6 d^2 (d = 54), so a class is multiplied out
 * only when N_r >= cross * d^2 (--cross, default 2); smaller
 * classes take bin/mom's own sequential steps on the sparse blocks, so in
 * that regime binmom costs what mom costs.
 *
 * Exact throughout: integer arithmetic, one division per entry at the end
 * of each class.  The final descent and the reports are bin/mom's
 * (mdecrease, perfindices), so the output is byte-identical to mom's.
 *
 * Prototype limits, each refused with a message rather than mis-solved:
 * mixed zero/non-zero think times (the lower replica levels of mom would
 * join the product state), and singular systems (use bin/mom, whose
 * util/singular.h fallback solves them exactly).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <gmpla.h>
#include "mom.h"
#include "lvtab.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* globals expected by mom's sources */
int AORSCTR, MULCTR, DIVCTR;
double AORSTIME, MULTIME, DIVTIME;
double t0, t1;
struct rusage ruse;
qnmodel* qnm;
mpq_t* MOM_CAPTURE = NULL;
bool INTERACTIVE, RANDGEN, CANON, ZSCALE, DEBUG, VERBOSE;
bool USE_LINBOX = false;
/* blocksolve calls this only when USE_LINBOX is set, which binmom never does */
int linbox_solve_dense(mpq_mat_t A, mpq_vec_t b, int N) { (void)A; (void)b; (void)N; return -1; }

static int    LEAF = 8;        /* steps multiplied directly at a leaf */
static double CROSS = 2.0;     /* binary splitting when N_r >= CROSS * d^2 */
static int    VERB = 0;
static int    NTHREADS = 1;

/* ---- dense integer matrices, row-major d x d ---- */
typedef struct { int d; mpz_t* a; } zmat;

static zmat zm_new(int d) { zmat m; int i; m.d = d; m.a = malloc((size_t)d*d*sizeof(mpz_t)); for (i=0;i<d*d;i++) mpz_init(m.a[i]); return m; }
static void zm_free(zmat m) { int i; for (i=0;i<m.d*m.d;i++) mpz_clear(m.a[i]); free(m.a); }

/* M(n) = n P' + Q' */
static zmat zm_step(const zmat* P, const zmat* Q, long n)
{
	zmat m = zm_new(P->d);
	int i;
	for (i=0;i<P->d*P->d;i++) { mpz_mul_si(m.a[i], P->a[i], n); mpz_add(m.a[i], m.a[i], Q->a[i]); }
	return m;
}

/* C = X Y */
static zmat zm_mul(const zmat* X, const zmat* Y)
{
	int d = X->d, i;
	zmat C = zm_new(d);
#ifdef _OPENMP
	#pragma omp parallel for num_threads(NTHREADS) schedule(dynamic,1) if(NTHREADS>1)
#endif
	for (i=0;i<d;i++) {
		int j, k;
		for (k=0;k<d;k++) {
			mpz_srcptr x = X->a[i*d+k];
			if (mpz_sgn(x) == 0) continue;
			for (j=0;j<d;j++)
				if (mpz_sgn(Y->a[k*d+j])) mpz_addmul(C.a[i*d+j], x, Y->a[k*d+j]);
		}
	}
	return C;
}

/* v <- X v */
static void zm_apply(const zmat* X, mpz_t* v)
{
	int d = X->d, i;
	mpz_t* w = malloc((size_t)d*sizeof(mpz_t));
	for (i=0;i<d;i++) mpz_init(w[i]);
#ifdef _OPENMP
	#pragma omp parallel for num_threads(NTHREADS) schedule(dynamic,1) if(NTHREADS>1)
#endif
	for (i=0;i<d;i++) {
		int k;
		for (k=0;k<d;k++)
			if (mpz_sgn(X->a[i*d+k]) && mpz_sgn(v[k])) mpz_addmul(w[i], X->a[i*d+k], v[k]);
	}
	for (i=0;i<d;i++) { mpz_swap(v[i], w[i]); mpz_clear(w[i]); }
	free(w);
}

/* M(b) ... M(a) */
static zmat zm_prod(const zmat* P, const zmat* Q, long a, long b)
{
	if (b - a + 1 <= LEAF) {
		zmat acc = zm_step(P, Q, a);
		long n;
		for (n=a+1;n<=b;n++) {
			zmat s = zm_step(P, Q, n), t = zm_mul(&s, &acc);
			zm_free(s); zm_free(acc); acc = t;
		}
		return acc;
	}
	long mid = a + (b - a) / 2;
	zmat lo = zm_prod(P, Q, a, mid), hi = zm_prod(P, Q, mid+1, b);
	zmat r = zm_mul(&hi, &lo);
	zm_free(lo); zm_free(hi);
	return r;
}

/* v <- M(b) ... M(a) v, sequentially or by binary splitting */
static void steps_apply(const zmat* P, const zmat* Q, long a, long b, mpz_t* v, int split)
{
	if (a > b) return;
	if (!split || b - a + 1 <= LEAF) {
		long n;
		for (n=a;n<=b;n++) { zmat s = zm_step(P, Q, n); zm_apply(&s, v); zm_free(s); }
		return;
	}
	long mid = a + (b - a) / 2;
	steps_apply(P, Q, a, mid, v, split);
	zmat hi = zm_prod(P, Q, mid+1, b);
	zm_apply(&hi, v);
	zm_free(hi);
}

/* ---- the class step as integer matrices ---- */

/* P' = D P and Q' = D Q from the class-r blocks of A; D returned.  0 on
 * success, -1 if a diagonal block is singular. */
static int class_matrices(LS* A, int cardGk, int cardG, zmat* Pz, zmat* Qz, mpz_t D)
{
	int d = cardGk + cardG, j, i;
	mpq_t* P = malloc((size_t)d*d*sizeof(mpq_t));
	mpq_t* Q = malloc((size_t)d*d*sizeof(mpq_t));
	for (i=0;i<d*d;i++) { mpq_init(P[i]); mpq_init(Q[i]); }
	mpq_vec_t e  = (mpq_vec_t) mpq_vec(d, 0, 1);
	mpq_vec_t b  = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	mpq_vec_t c  = (mpq_vec_t) mpq_vec(cardG, 0, 1);
	mpq_vec_t b2 = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	int rc = 0;
	for (j=0;j<d && rc==0;j++) {
		mpq_set_ui(e[j], 1, 1);
		/* column j of P: [A11^{-1} B1r e_j; 0] */
		mpq_mspvecmul(b, A->B1r, e);
		mpq_vec_t x = blocksolve(A, b);
		if (!x) { rc = -1; break; }
		for (i=0;i<cardGk;i++) mpq_set(P[i*d+j], x[i]);
		for (i=0;i<cardGk;i++) mpq_clear(x[i]);
		free(x);
		/* column j of Q: [-A11^{-1} A12 B2r e_j; B2r e_j] */
		mpq_mspvecmul(c, A->B2r, e);
		for (i=0;i<cardG;i++) mpq_set(Q[(cardGk+i)*d+j], c[i]);
		mpq_mspvecmul(b2, A->A12, c);
		x = blocksolve(A, b2);
		if (!x) { rc = -1; break; }
		for (i=0;i<cardGk;i++) mpq_neg(Q[i*d+j], x[i]);
		for (i=0;i<cardGk;i++) mpq_clear(x[i]);
		free(x);
		mpq_set_ui(e[j], 0, 1);
	}
	if (rc == 0) {
		mpz_set_ui(D, 1);
		for (i=0;i<d*d;i++) { mpz_lcm(D, D, mpq_denref(P[i])); mpz_lcm(D, D, mpq_denref(Q[i])); }
		*Pz = zm_new(d); *Qz = zm_new(d);
		for (i=0;i<d*d;i++) {
			mpz_divexact(Pz->a[i], D, mpq_denref(P[i])); mpz_mul(Pz->a[i], Pz->a[i], mpq_numref(P[i]));
			mpz_divexact(Qz->a[i], D, mpq_denref(Q[i])); mpz_mul(Qz->a[i], Qz->a[i], mpq_numref(Q[i]));
		}
	}
	for (i=0;i<d*d;i++) { mpq_clear(P[i]); mpq_clear(Q[i]); }
	free(P); free(Q);
	for (i=0;i<=d;i++) mpq_clear(e[i]);
	for (i=0;i<=cardGk;i++) { mpq_clear(b[i]); mpq_clear(b2[i]); }
	for (i=0;i<=cardG;i++) mpq_clear(c[i]);
	free(e); free(b); free(b2); free(c);
	return rc;
}

/* g (rational, d) <- g(N_r) and gr <- g(N_r - 1) for class r */
static void class_run(const zmat* P, const zmat* Q, const mpz_t D, long N, int split, mpq_vec_t g, mpq_vec_t gr)
{
	int d = P->d, i;
	mpz_t den0, f;
	mpz_init_set_ui(den0, 1); mpz_init(f);
	for (i=0;i<d;i++) mpz_lcm(den0, den0, mpq_denref(g[i]));
	mpz_t* v = malloc((size_t)d*sizeof(mpz_t));
	for (i=0;i<d;i++) { mpz_init(v[i]); mpz_divexact(v[i], den0, mpq_denref(g[i])); mpz_mul(v[i], v[i], mpq_numref(g[i])); }

	steps_apply(P, Q, 1, N-1, v, split);      /* g(N-1), unscaled */
	/* gr = v / (den0 D^(N-1) (N-1)!) */
	mpz_pow_ui(f, D, (unsigned long)(N-1));
	{ mpz_t fa; mpz_init(fa); mpz_fac_ui(fa, (unsigned long)(N-1)); mpz_mul(f, f, fa); mpz_clear(fa); }
	mpz_mul(f, f, den0);
	for (i=0;i<d;i++) { mpz_set(mpq_numref(gr[i]), v[i]); mpz_set(mpq_denref(gr[i]), f); mpq_canonicalize(gr[i]); }

	steps_apply(P, Q, N, N, v, 0);           /* last step */
	mpz_mul(f, f, D); mpz_mul_ui(f, f, (unsigned long)N);
	for (i=0;i<d;i++) { mpz_set(mpq_numref(g[i]), v[i]); mpz_set(mpq_denref(g[i]), f); mpq_canonicalize(g[i]); }

	for (i=0;i<d;i++) mpz_clear(v[i]);
	free(v); mpz_clear(den0); mpz_clear(f);
}

/* bin/mom's sequential steps on the sparse blocks, for classes too short
 * to repay forming the dense d x d step matrices (N_r < cross * d).  0, or
 * -1 if a block solve fails. */
static int class_seq(LS* A, int cardGk, int cardG, long N, mpq_vec_t g, mpq_vec_t gr)
{
	int t;
	long n;
	mpq_vec_t G  = (mpq_vec_t) mpq_vec(cardG, 0, 1);
	mpq_vec_t b  = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	mpq_vec_t b2 = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	mpq_t ninv; mpq_init(ninv);
	int rc = 0;
	for (n=1;n<=N;n++) {
		if (n == N) mpq_vecdup(gr, g, cardGk+cardG);
		mpq_mspvecmul(G, A->B2r, g);
		mpq_set_si(ninv, 1, (unsigned long) n);
		for (t=0;t<cardG;t++) mpq_mul(G[t], ninv, G[t]);
		mpq_mspvecmul(b, A->B1r, g);
		mpq_mspvecmul(b2, A->A12, G);
		for (t=0;t<cardGk;t++) mpq_sub(b[t], b[t], b2[t]);
		mpq_vec_t Gk = blocksolve(A, b);
		if (!Gk) { rc = -1; break; }
		for (t=0;t<cardGk;t++) mpq_set(g[t], Gk[t]);
		for (t=0;t<cardG;t++)  mpq_set(g[cardGk+t], G[t]);
		for (t=0;t<=cardGk;t++) mpq_clear(Gk[t]);
		free(Gk);
	}
	mpq_clear(ninv);
	for (t=0;t<=cardG;t++) mpq_clear(G[t]);
	for (t=0;t<=cardGk;t++) { mpq_clear(b[t]); mpq_clear(b2[t]); }
	free(G); free(b); free(b2);
	return rc;
}

static void freels(LS* A, int m)
{
	int t = 0, h, k;
	if (!A) return;
	for (h=1;h<=A->H;h++) for (k=0;k<nck(m,h);k++) { free(A->C[t].diag); free(A->C[t++].nondiag->coeff); }
	free(A->A12); free(A->B1r); free(A->B2r); free(A);
}

static void vfree(mpq_vec_t v, int n) { int i; if (!v) return; for (i=0;i<=n;i++) mpq_clear(v[i]); free(v); }

static int solve(qnmodel* q, bool vb, bool lo, bool ne, bool ng, bool tp, bool ql, bool db, mpz_t sf)
{
	int m = q->M, R = q->R, r, i, s, t, cardG = 0, cardGk = 0, dlen = 0;
	int* n = (int*) int_vec(R, 0);
	LS* A = NULL;
	mpq_vec_t g = NULL, gr = NULL;
	nckinit(m+R-1, R);

	for (r=1;r<=R;r++) {
		double tc = CPUTIME;
		cardGk = nck(m+r-1, r)*r;
		cardG  = nck(m+r-2, r-1)*r;
		int d = cardGk + cardG;
		freels(A, m);
		A = setupls(q->L, n, q->Z, q->mi, m, r);
		if (!A) { fprintf(stderr, "binmom: singular system at class %d; use bin/mom (exact fallback)\n", r); return 1; }
		if (r == 1) {
			g  = (mpq_vec_t) mpq_vec(m+1, 1, 1);
			gr = (mpq_vec_t) mpq_vec(m+1, 0, 1);
		} else {
			/* class transition, as in bin/mom */
			mpq_vec_t G  = (mpq_vec_t) mpq_vec(cardG, 0, 1);
			mpq_vec_t b  = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
			for (i=0;i<nck(m+r-2, r-1);i++) {
				for (s=0;s<=r-2;s++) mpq_set(G[i*r+s], g[i*(r-1)+s]);
				mpq_set(G[i*r+r-1], gr[i*(r-1)]);
			}
			mpq_mspvecmul(b, A->A12, G);
			for (t=0;t<cardGk;t++) mpq_neg(b[t], b[t]);
			mpq_vec_t Gk = blocksolve(A, b);
			if (!Gk) { fprintf(stderr, "binmom: singular system at class %d; use bin/mom (exact fallback)\n", r); return 1; }
			vfree(g, dlen); vfree(gr, dlen);
			g  = (mpq_vec_t) mpq_vec(d, 0, 1);
			gr = (mpq_vec_t) mpq_vec(d, 0, 1);
			for (t=0;t<cardGk;t++) mpq_set(g[t], Gk[t]);
			for (t=0;t<cardG;t++)  mpq_set(g[cardGk+t], G[t]);
			vfree(Gk, cardGk); vfree(G, cardG); vfree(b, cardGk);
		}
		dlen = d;
		if (q->N[r-1] > 0 && (double) q->N[r-1] < CROSS * (double) d * d) {
			double ts = CPUTIME;
			if (class_seq(A, cardGk, cardG, q->N[r-1], g, gr) < 0) {
				fprintf(stderr, "binmom: singular system at class %d; use bin/mom (exact fallback)\n", r);
				return 1;
			}
			if (VERB)
				fprintf(stderr, "binmom: class %d  d=%d  N_r=%d  sequential (mom steps)  steps %.3fs\n",
				        r, d, q->N[r-1], CPUTIME - ts);
		} else if (q->N[r-1] > 0) {
			zmat P, Q; mpz_t D; mpz_init(D);
			if (class_matrices(A, cardGk, cardG, &P, &Q, D) < 0) {
				fprintf(stderr, "binmom: singular system at class %d; use bin/mom (exact fallback)\n", r);
				return 1;
			}
			double ts = CPUTIME;
			class_run(&P, &Q, D, q->N[r-1], 1, g, gr);
			if (VERB)
				fprintf(stderr, "binmom: class %d  d=%d  N_r=%d  %s  D %lu bits  setup %.3fs  steps %.3fs\n",
				        r, d, q->N[r-1], "binary splitting",
				        (unsigned long) mpz_sizeinbase(D, 2), ts - tc, CPUTIME - ts);
			zm_free(P); zm_free(Q); mpz_clear(D);
		}
		n[r-1] = q->N[r-1];
	}
	{
		mpq_vec_t G = (mpq_vec_t) mpq_vec(1, 0, 1), Gk = (mpq_vec_t) mpq_vec(1, 0, 1);
		mdecrease(q, G, Gk, g, gr, vb, lo, ne, ng, tp, ql, db, false, sf);
	}
	freels(A, m);
	free(n);
	return 0;
}

static void usage(const char* p)
{
	printf("USAGE: %s [-l|-e|-g|-t|-q|-d] [-v] [-j n] [--cross c] [--leaf k] [--seq] model.qn\n", p);
	printf("  binmom: MoM whose population steps are multiplied together by binary\n");
	printf("  splitting (exact, integer arithmetic); output identical to bin/mom.\n");
	printf("  -l -e -g -t -q -d : as bin/mom\n");
	printf("  -v         : per-class basis size, route and timing on stderr\n");
	printf("  -j n       : threads for the matrix products (OpenMP)\n");
	printf("  --cross c  : binary splitting for a class when N_r >= c * d^2 (default 2)\n");
	printf("  --leaf k   : steps multiplied directly at a leaf (default 8)\n");
	printf("  --seq      : never split (bin/mom's sequential steps throughout)\n");
	printf("  Not handled (refused): mixed zero/non-zero think times, singular systems.\n");
}

int main(int argc, char** argv)
{
	bool lo=false, ne=false, ng=false, tp=false, ql=false, db=false;
	char* file = NULL;
	int i;
	t0 = CPUTIME;
	for (i=1;i<argc;i++) {
		if      (!strcmp(argv[i],"-l")||!strcmp(argv[i],"--log"))   lo = true;
		else if (!strcmp(argv[i],"-e")||!strcmp(argv[i],"--ex"))    ne = true;
		else if (!strcmp(argv[i],"-g")||!strcmp(argv[i],"--nc"))    ng = true;
		else if (!strcmp(argv[i],"-t")||!strcmp(argv[i],"--tput"))  tp = true;
		else if (!strcmp(argv[i],"-q")||!strcmp(argv[i],"--qlen"))  ql = true;
		else if (!strcmp(argv[i],"-d")||!strcmp(argv[i],"--debug")) db = true;
		else if (!strcmp(argv[i],"-v")) VERB = 1;
		else if (!strcmp(argv[i],"-j") && i+1<argc) NTHREADS = atoi(argv[++i]);
		else if (!strcmp(argv[i],"--cross") && i+1<argc) CROSS = atof(argv[++i]);
		else if (!strcmp(argv[i],"--leaf") && i+1<argc) { LEAF = atoi(argv[++i]); if (LEAF < 1) LEAF = 1; }
		else if (!strcmp(argv[i],"--seq")) CROSS = 1e300;
		else if (!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")) { usage(argv[0]); return 0; }
		else if (argv[i][0] != '-') file = argv[i];
	}
	if (!file) { usage(argv[0]); return -1; }
	qnm = readmodel(file);
	qn_require_closed(qnm, "binmom", 0);
	if (mom_mixed_think(qnm)) {
		fprintf(stderr, "binmom: mixed zero/non-zero think times are not handled by this prototype; use bin/mom\n");
		return 1;
	}
	mpz_t sf; mpz_init_set_ui(sf, 1);
	bool report = !lo && !ne && !ng && !tp && !ql;
	if (report) printmodel(qnm);
	int rc = solve(qnm, false, lo, ne, ng, tp, ql, db, sf);
	if (rc == 0 && report) { t1 = CPUTIME; printf("Elapsed time (binmom): %.6f s\n", t1 - t0); }
	mpz_clear(sf);
	return rc;
}
