#ifndef UTIL 
#define UTIL

#ifndef __cplusplus
  #ifndef _STDBOOL_H
    #define bool unsigned int
    #define true 1
    #define false 0
  #endif
#endif
#define MIN(a,b) ((a<b) ? a : b)
#define MAX(a,b) ((a>b) ? a : b)
#define MAXNCKTABLE 100 /* maximum allowed factorial */

typedef struct 
{
	int M; /* number of queues */
	int R; /* number of classes */
	int *N; /* job populations */
	int *Z; /* think times */
	int **L; /* service demands */
	int *mi; /* multiplicities */
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
void printmodel(qnmodel* qn);

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
