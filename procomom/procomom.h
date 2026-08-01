#ifndef PROCOMOM_H
#define PROCOMOM_H

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <gmpla.h>
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

/* Combinations with repetition (independent from comom) */
typedef struct
{
	int n;       /* number of elements (R classes) */
	int k;       /* choose from k types (M stations) */
	int card;    /* nchoosek(n+k-1, k) */
	int **combs; /* matrix of combinations */
	int stride;  /* components per shift: M normally (base + stations 1..M-1),
	                M+1 when the reference station is replicated and its own
	                replica term needs a column of its own.  See phash.c. */
} combsrep;

/* P matrices: 4 dense rectangular matrices [numRows x basisSize] */
typedef struct
{
	mpq_mat_t A;
	mpq_mat_t B;
	mpq_mat_t DC;
	mpq_mat_t DD;
	int numRows;
	int basisSize;
} PMatrices;

extern qnmodel* qnm;
extern double t0, t1;
extern struct rusage ruse;

/* Function prototypes */
int phash(combsrep* Dn, int* comb, int i);
PMatrices* genpmatrix(combsrep* Dn, qnmodel* qnm, int* Ncur, int r);
void pexact(mpq_vec_t pk, combsrep* Dn, int M);

#endif
