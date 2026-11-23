/*
 * binsplit.h - exact products of the population steps of one class by
 * binary splitting, shared by binmom and bincomom (header-only).
 *
 * Both recursions advance the basis g of a class by
 *
 *   n g(n) = (n P + Q) g(n-1),   P, Q fixed within the class,
 *
 * so with P' = D P, Q' = D Q integer (D the common denominator)
 *
 *   g(b) = (b P' + Q') ... (a P' + Q') g(a-1) / (D^(b-a+1) b!/(a-1)!).
 *
 * The product is formed in a balanced tree, so the large multiplications
 * are between integers of equal size and GMP's subquadratic multiplication
 * applies.  BS_LEAF steps are multiplied directly at a leaf; BS_NTHREADS
 * threads share the rows of each matrix product.
 */
#ifndef BINSPLIT_H
#define BINSPLIT_H
#include <stdlib.h>
#include <gmp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int BS_LEAF = 8;
static int BS_NTHREADS = 1;

/* dense integer matrix, row-major d x d */
typedef struct { int d; mpz_t* a; } zmat;

static inline zmat zm_new(int d) { zmat m; int i; m.d = d; m.a = malloc((size_t)d*d*sizeof(mpz_t)); for (i=0;i<d*d;i++) mpz_init(m.a[i]); return m; }
static inline void zm_free(zmat m) { int i; for (i=0;i<m.d*m.d;i++) mpz_clear(m.a[i]); free(m.a); }

/* M(n) = n P' + Q' */
static inline zmat zm_step(const zmat* P, const zmat* Q, long n)
{
	zmat m = zm_new(P->d);
	int i;
	for (i=0;i<P->d*P->d;i++) { mpz_mul_si(m.a[i], P->a[i], n); mpz_add(m.a[i], m.a[i], Q->a[i]); }
	return m;
}

/* C = X Y */
static inline zmat zm_mul(const zmat* X, const zmat* Y)
{
	int d = X->d, i;
	zmat C = zm_new(d);
#ifdef _OPENMP
	#pragma omp parallel for num_threads(BS_NTHREADS) schedule(dynamic,1) if(BS_NTHREADS>1)
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
static inline void zm_apply(const zmat* X, mpz_t* v)
{
	int d = X->d, i;
	mpz_t* w = malloc((size_t)d*sizeof(mpz_t));
	for (i=0;i<d;i++) mpz_init(w[i]);
#ifdef _OPENMP
	#pragma omp parallel for num_threads(BS_NTHREADS) schedule(dynamic,1) if(BS_NTHREADS>1)
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
static inline zmat zm_prod(const zmat* P, const zmat* Q, long a, long b)
{
	if (b - a + 1 <= BS_LEAF) {
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

/* v <- M(b) ... M(a) v, step by step or by binary splitting */
static inline void bs_steps_apply(const zmat* P, const zmat* Q, long a, long b, mpz_t* v, int split)
{
	if (a > b) return;
	if (!split || b - a + 1 <= BS_LEAF) {
		long n;
		for (n=a;n<=b;n++) { zmat s = zm_step(P, Q, n); zm_apply(&s, v); zm_free(s); }
		return;
	}
	long mid = a + (b - a) / 2;
	bs_steps_apply(P, Q, a, mid, v, split);
	zmat hi = zm_prod(P, Q, mid+1, b);
	zm_apply(&hi, v);
	zm_free(hi);
}

/* P' = D P, Q' = D Q from rational d x d row-major P, Q; D set to the lcm
 * of all denominators. */
static inline void bs_integerize(const mpq_t* P, const mpq_t* Q, int d, zmat* Pz, zmat* Qz, mpz_t D)
{
	int i;
	mpz_set_ui(D, 1);
	for (i=0;i<d*d;i++) { mpz_lcm(D, D, mpq_denref(P[i])); mpz_lcm(D, D, mpq_denref(Q[i])); }
	*Pz = zm_new(d); *Qz = zm_new(d);
	for (i=0;i<d*d;i++) {
		mpz_divexact(Pz->a[i], D, mpq_denref(P[i])); mpz_mul(Pz->a[i], Pz->a[i], mpq_numref(P[i]));
		mpz_divexact(Qz->a[i], D, mpq_denref(Q[i])); mpz_mul(Qz->a[i], Qz->a[i], mpq_numref(Q[i]));
	}
}

/* Integer form of a rational vector: v[i] = g[i] * den, den = lcm of the
 * denominators of g.  v must hold d initialized mpz_t. */
static inline void bs_vec_to_z(const mpq_t* g, int d, mpz_t* v, mpz_t den)
{
	int i;
	mpz_set_ui(den, 1);
	for (i=0;i<d;i++) mpz_lcm(den, den, mpq_denref(g[i]));
	for (i=0;i<d;i++) { mpz_divexact(v[i], den, mpq_denref(g[i])); mpz_mul(v[i], v[i], mpq_numref(g[i])); }
}

/* g[i] = v[i] / den, canonical */
static inline void bs_vec_from_z(mpq_t* g, int d, const mpz_t* v, const mpz_t den)
{
	int i;
	for (i=0;i<d;i++) { mpz_set(mpq_numref(g[i]), v[i]); mpz_set(mpq_denref(g[i]), den); mpq_canonicalize(g[i]); }
}

/* Scale den by D^(b-a+1) b!/(a-1)!: the divisor the steps a..b add. */
static inline void bs_scale_den(mpz_t den, const mpz_t D, long a, long b)
{
	mpz_t f, h;
	if (a > b) return;
	mpz_init(f); mpz_init(h);
	mpz_pow_ui(f, D, (unsigned long)(b - a + 1));
	mpz_mul(den, den, f);
	mpz_fac_ui(f, (unsigned long) b);
	mpz_fac_ui(h, (unsigned long)(a - 1));
	mpz_divexact(f, f, h);
	mpz_mul(den, den, f);
	mpz_clear(f); mpz_clear(h);
}

#endif
