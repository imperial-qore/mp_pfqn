#ifndef UTIL 
#define UTIL

#include <gmp.h>

/* Use the standard bool unconditionally.
 *
 * The previous fallback defined bool as unsigned int whenever <stdbool.h>
 * had not already been included in the translation unit.  Since some
 * sources include <stdbool.h> and others do not, the same global bool was
 * a 1-byte _Bool in one object and a 4-byte unsigned int in another.
 * Reading such a global through the wider type picks up three adjacent
 * bytes of unrelated data, so its value varied between runs: in mom this
 * made USE_LINBOX read as true at random in setupls.c/blocksolve.c,
 * silently taking a different solve path from the one that built the
 * factors and producing wrong results (including negative queue lengths)
 * in a few percent of runs. */
#ifndef __cplusplus
  #include <stdbool.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b) ((a<b) ? a : b)
#define MAX(a,b) ((a>b) ? a : b)
#define MAXNCKTABLE 100 /* maximum allowed factorial */

typedef struct
{
	int M; /* number of queues */
	int R; /* number of classes */
	int *N; /* job populations, N[r]=-1 means open class */
	mpz_t *Z; /* think times */
	mpz_t **L; /* service demands */
	int *mi; /* multiplicities */
	int hasOpen; /* 1 if LAMBDA section present (open classes) */
	mpq_t *lambda; /* arrival rates [R], NULL if no open classes */
	int isLD; /* 1 if MU section present */
	int Nt; /* total closed population */
	mpq_t **mu; /* load-dependent rates [M][Nt], NULL if not LD */
} qnmodel;

extern int nckmaxn; /* nchoosek table */
extern int nckmaxk; /* nchoosek table */

/* general math */
long int factorial(long int N);
int ** multichoose(int n,int k);
long long int nck(int N, int K);
void nckinit(int nmax, int kmax);
int nnz(int* v, int n);
int nnzposcmp(int*i1, int*i2, int n);
int randi(int min,int max);
int** sortbynnzpos(int** I, int m,int n);
int sum(int* v, int n);

/* queueing network model */
qnmodel* readmodel(char* filename);
void freemodel(qnmodel* qn);
void printmodel(qnmodel* qn);
void printmodel_with_perturbation(qnmodel* qn, int perturbation_digit, mpz_t scale_factor, int perturbation_seed, mpz_t** original_L, mpz_t* original_Z);
int* getclosedclasses(qnmodel* qn, int* numClosed);
int* getopenclasses(qnmodel* qn, int* numOpen);

/* queueing network population */
int* initpop(int R);
void printpop(int*n, int R);
int* resetpop(int*,int R);
int nextpop(int* n, int* N, int R);
int nextpopinteractive(int* n, int* N, int R);
int popindex(int* n,int R,int* planesizes);
bool issingleclass(int* n, int R);
int* getplanesizes(int *N, int R);

/* log of a positive rational without converting it to a double first:
 * mpq_get_d overflows to inf once the value exceeds ~1e308, which a
 * normalizing constant does at modest populations.  mpz_get_d_2exp splits
 * off the binary exponent of numerator and denominator separately. */
static inline double mpq_log(const mpq_t q)
{
	/* When q fits in a double, log(mpq_get_d(q)) rounds once and gives
	 * the same last digit as the other solvers' output; the two-term
	 * split below rounds twice, so use it only when q over/underflows. */
	double d = mpq_get_d(q);
	if (isfinite(d) && d > 0.0 && isnormal(d)) return log(d);
	signed long en, ed;
	double dn = mpz_get_d_2exp(&en, mpq_numref(q));
	double dd = mpz_get_d_2exp(&ed, mpq_denref(q));
	return log(dn) - log(dd) + (double)(en - ed) * log(2.0);
}

/* Same for an mpf_t (log(mpf_get_d(x)) overflows the same way). */
static inline double mpf_log(const mpf_t x)
{
	double v = mpf_get_d(x);
	if (isfinite(v) && v > 0.0 && isnormal(v)) return log(v); /* single rounding */
	signed long e;
	double d = mpf_get_d_2exp(&e, x);
	return log(d) + (double)e * log(2.0);
}

/* Solvers that only handle closed product-form models call this right
 * after readmodel: an open class (N[r]<0) gives a population lattice of
 * size zero and a load-dependent MU section is silently ignored, so
 * either would otherwise print a wrong answer (or crash) with exit 0. */
static inline void qn_require_closed(const qnmodel* qn, const char* solver, int allow_ld)
{
	int r;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] < 0) {
			fprintf(stderr, "%s: class %d is open (N=%d); this solver supports closed models only\n", solver, r+1, qn->N[r]);
			exit(1);
		}
	}
	if (qn->hasOpen) {
		fprintf(stderr, "%s: LAMBDA section present; this solver supports closed models only\n", solver);
		exit(1);
	}
	if (qn->isLD && !allow_ld) {
		fprintf(stderr, "%s: MU (load-dependent) section present; this solver does not support load-dependent stations\n", solver);
		exit(1);
	}
}

#endif
