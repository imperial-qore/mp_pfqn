#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h" 

long long int ncktable[MAXNCKTABLE][MAXNCKTABLE]; /* nchoosek table */

long long int nck(int N, int K)
{
	// Use symmetry: nck(n,k) = nck(n,n-k)
	// This avoids overflow for large k values
	if (K > N - K) {
		K = N - K;
	}
	
	if (N>nckmaxn || K>nckmaxk || ncktable[N][K]==0 ) 
	{
	long int n=N;
	long int k=K;
	long int i;
	long long int result;
	if(k==0) 
	{
		ncktable[N][K]=1;
		return 1;
	}
	result = n;
	i=2;
	n--;k--;
	while(k>0)
	{
		result = (long long int) (result*n)/i;
		i++;
		n--;k--;
	}
	ncktable[N][K]=result;
	return result;
	}
	else
	{
	return ncktable[N][K];
	}
}
