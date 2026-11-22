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

#endif
