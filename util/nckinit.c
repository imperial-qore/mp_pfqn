#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include <time.h>
#include "util.h"

/* Define the global variables here */
int nckmaxn;
int nckmaxk;

void nckinit(int nmax, int kmax)
{
//	nck=(int**)int_matzeros(nmax+1,kmax+1);
	nckmaxn=nmax;
	nckmaxk=kmax;
//	int i,j;
//	for (i=0;i<=nmax;i++)
//		for (j=0;j<=kmax;j++)
//			nck[i][j]=0;
}
