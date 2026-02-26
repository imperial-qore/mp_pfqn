#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "gmpla.h"

int* mpq_ludcmp(mpq_mat_t A, int N)
{
 int i,j,k;
 int imax=0;
 int* indx=(int*)calloc(N+1,sizeof(int));
 mpq_t max; mpq_init(max);
 mpq_t sum; mpq_init(sum);
 mpq_t t; mpq_init(t);
 mpq_t w; mpq_init(w);
 mpq_vec_t vv=(mpq_vec_t)mpq_vec(N+1,0,1);

 for(i=1;i<=N;i++)
	 mpq_init(vv[i-1]);
 mpq_set_d(t,1); // serve nel loop
 for(i=1;i<=N;i++)
 {
 	mpq_set_d(max,0);
	for (j=1; j<=N; j++)
	{
		mpq_abs(w,A[i-1][j-1]);
		if (mpq_cmp(w,max)>0)
			mpq_set(max,w);
	}
	if (mpq_sgn(max)==0) /* singular matrix */
	{
		// Free allocated memory before returning
		for(i=0;i<N;i++)
			mpq_clear(vv[i]);
		free(vv);
		free(indx);
		mpq_clear(max);
		mpq_clear(sum);
		mpq_clear(t);
		mpq_clear(w);
		return NULL;
	}
	mpq_div(vv[i-1],t,max);
 }
 for(j=1;j<=N;j++)
 {
	int i;
 	for(i=1;i<=j-1;i++)
 	{
		mpq_set(sum, A[i-1][j-1]);
 		for(k=1;k<=i-1;k++)
		{
		  if (mpq_sgn(A[i-1][k-1])!=0 && mpq_sgn(A[k-1][j-1])!=0) /* singular matrix */
		  {
			mpq_mul(t,A[i-1][k-1],A[k-1][j-1]);
			mpq_sub(sum,sum,t);
		  }
		}
		mpq_set(A[i-1][j-1],sum);
	}
	mpq_set_d(max,0);
 	for( i=j;i<=N;i++)
	{
		mpq_set(sum, A[i-1][j-1]);
		for(k=1;k<=j-1;k++)
		{
		  if (mpq_sgn(A[i-1][k-1])!=0 && mpq_sgn(A[k-1][j-1])!=0) /* singular matrix */
		  {
		  mpq_mul(t,A[i-1][k-1],A[k-1][j-1]);
		  mpq_sub(sum,sum,t);
		  }
		}
		mpq_set(A[i-1][j-1],sum);
		mpq_abs(w,sum); 
		mpq_mul(w,w,vv[i-1]);
		if (mpq_cmp(w,max)>=0)
		{
			mpq_set(max,w);
			imax=i;
		}
	}
	if (j!=imax)
	{
 		for(k=1;k<=N;k++)
		{
			mpq_set(w,A[imax-1][k-1]);
			mpq_set(A[imax-1][k-1],A[j-1][k-1]);
			mpq_set(A[j-1][k-1],w);
		}
		mpq_set(vv[imax-1],vv[j-1]);
	}
	indx[j-1]=imax;
	if (j!=N)
	{
		mpq_set_d(t,1.0);
		if (mpq_sgn(A[j-1][j-1])==0) /* singular matrix */
		{
			// Free allocated memory before returning
			for(i=0;i<N;i++)
				mpq_clear(vv[i]);
			free(vv);
			free(indx);
			mpq_clear(max);
			mpq_clear(sum);
			mpq_clear(t);
			mpq_clear(w);
			return NULL;
		}

		mpq_div(w,t,A[j-1][j-1]);
 		for(i=j+1;i<=N;i++)
		{
		    if (mpq_sgn(A[i-1][j-1])!=0 && mpq_sgn(w)!=0)
		    {
			mpq_mul(A[i-1][j-1],A[i-1][j-1],w);
		    }
		}
	}
 }
 /* Check last pivot (the main loop skips this check when j==N) */
 if (mpq_sgn(A[N-1][N-1])==0) /* singular matrix */
 {
	for(i=0;i<N;i++)
		mpq_clear(vv[i]);
	free(vv);
	free(indx);
	mpq_clear(max);
	mpq_clear(sum);
	mpq_clear(t);
	mpq_clear(w);
	return NULL;
 }
 mpq_clear(max);
 mpq_clear(sum);
 mpq_clear(t);
 mpq_clear(w);
 for(i=0;i<N;i++)
	mpq_clear(vv[i]);
 free(vv);
 return indx;
}
