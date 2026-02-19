#include <stdbool.h>

int* initpop(int R);
void printpop(int*n, int R);
int* resetpop(int*,int R);
int nextpop(int* n, int* N, int R);
int nextpopinteractive(int* n, int* N, int R);
int popindex(int* n,int R,int* planesizes);
bool issingleclass(int* n, int R);
int* getplanesizes(int *N, int R);
void convolution_multi_bignum(mpf_t *g, mpf_t *X);
