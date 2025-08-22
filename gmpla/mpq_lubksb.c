#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gmpla.h"

int mpq_lubksb(mpq_mat_t A, mpq_vec_t b,int N,int* indx)
{
	int i,j;
	int l=0; 
	mpq_t t; mpq_init(t);
	mpq_t sum; mpq_init(sum);
	/* solve Ly=b*/
	for (i =1; i<=N; i++) 
	{
	 int p = indx[i-1]; 
    	 mpq_set(sum,b[p-1]);
         mpq_set(b[p-1],b[i-1]); 
         if (l != 0)     
	 {
         for (j=l; j<=i-1;j++)
	 {
	 if (mpq_sgn(A[i-1][j-1])!=0)  
	 {
		mpq_mul(t,A[i-1][j-1],b[j-1]);
             	mpq_sub(sum,sum,t);
	 }
	 }
	 }
         else if (mpq_sgn(sum)!= 0)
	 {
            l = i;
     	 }
         mpq_set(b[i-1],sum);
	}
	/* solve Ux=y */
	for (i = N; i>=1; i--)
	{
  	 mpq_set(sum, b[i-1]);
	  for (j=i+1;j<=N;j++)
	  {
		if (mpq_sgn(A[i-1][j-1])!=0)  
		{
		mpq_mul(t,A[i-1][j-1],b[j-1]);

             	mpq_sub(sum,sum,t); /* BOTTLENECK OPERATION */
		}
	  }
    	  mpq_set(b[i-1],sum);
  	  if (mpq_sgn(A[i-1][i-1])==0) 
	  {
//		mpq_matprint(A,N,N);
		mpq_clear(t);
		mpq_clear(sum);
		return -1;  // Return error code for singular matrix
	  }
	  mpq_div(b[i-1],b[i-1],A[i-1][i-1]); 
	}
	mpq_clear(t);
	mpq_clear(sum);
	return 0;  // Success
}

/* test

int main()
{
	int i,j,n=4;
	mpq_t** A=(mpq_t**)mpq_matzeros(n,n);
	mpq_t*  b=(mpq_t* )mpq_veczeros(n);

	srand((unsigned)time(NULL));
	
	for (i=1;i<=n;i++)
	   for (j=1;j<=n;j++)
	   {
	      mpq_set_ui(A[i-1][j-1],(int)(rand()/1e6),i*j);
	      mpq_set_ui(b[i-1],(int)(rand()/1e6),i*j);
	   }

	mpq_matprint(A,n,n);
	mpq_vecprint(b,n);
	int* indices = mpq_ludcmp(A,n);
	if (indices==NULL) printf("ERROR: Singular matrix\n");
	mpq_matprint(A,n,n);
	mpq_lubksb(A,b,n,indices);
	mpq_vecprint(b,n);
	return 1;
}
*/
