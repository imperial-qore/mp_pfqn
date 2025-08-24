#include <stdio.h>
#include <stdlib.h>
#include <math.h>		
#include <gmp.h>
#include "mom.h"
#include "popcycle.h"
//#define DEBUG
void convolution_multi_bignum(mpf_t *g, mpf_t *X)
{
	int M=qnm->M;
	int R=qnm->R;
	int r,t,m;
	int*N=qnm->N;
	int *n=(int*)initpop(R);
  	int *planesizes=(int*) getplanesizes(N,R); // planesize[k]=\prod_{j=R..k} (Nj+1)
	
	// Check if any think times are non-zero
	int hasZ = 0;
	for (r=1; r<=R; r++) {
		if (qnm->Z[r-1] > 0) {
			hasZ = 1;
			break;
		}
	}
	
	// If think times present, need M+1 stations (including delay server)
	int Mz = hasZ ? M + 1 : M;
	
	// init: the storage space
	mpf_t** G=(mpf_t**) calloc(planesizes[0],sizeof(mpf_t*)); // the recursion is performed on N, so we keep the M queues of population N-1
	//mpf_set_default_prec(1024);
	//printf("convolution_multi_bignum: current precision is %ld\n",mpf_get_default_prec());
	
	for (t=1;t<=planesizes[0];t++)
	{
  		G[t-1]=(mpf_t*) calloc(Mz,sizeof(mpf_t));
		for (m=1;m<=Mz;m++)
		{
  			mpf_init(G[t-1][m-1]);
			mpf_set_d(G[t-1][m-1], (double) 0.0);
		}
	}
	for (m=1;m<=Mz;m++)
	{
		mpf_set_d(G[0][m-1], (double) 1.0);
	}
 	// end 
	int curindex=-1;
	do	
	{
	curindex++; // current population index
	
	// If think times present, first compute G[0,n] based on Z (delay server)
	if (hasZ) {
		mpf_set_d(G[curindex][0], (double) 1.0);
		for (r=1;r<=R;r++)
		{
			if (n[r-1]!=0 && qnm->Z[r-1] > 0)
			{
				// G[0,n] = prod(Z[r]^n[r] / n[r]!)
				mpf_t zr; mpf_init(zr);
				mpf_set_d(zr, (double) qnm->Z[r-1]);
				mpf_t zr_pow; mpf_init(zr_pow);
				mpf_pow_ui(zr_pow, zr, n[r-1]);
				
				// Divide by n[r]!
				mpf_t fact; mpf_init(fact);
				mpf_set_ui(fact, 1);
				for (int i = 1; i <= n[r-1]; i++) {
					mpf_mul_ui(fact, fact, i);
				}
				mpf_div(zr_pow, zr_pow, fact);
				
				mpf_mul(G[curindex][0], G[curindex][0], zr_pow);
				mpf_clear(zr);
				mpf_clear(zr_pow);
				mpf_clear(fact);
			}
		}
	}
	
	for (m=1;m<=M;m++)
	{
		int m_idx = hasZ ? m : m-1; // Adjust index when delay server is present
		
		if (m>1) // add norm const of the complementary network
		{
			mpf_set(G[curindex][m_idx],G[curindex][m_idx-1]);
		}
		else if (hasZ) // m==1 with think times, copy from delay server
		{
			mpf_set(G[curindex][m_idx],G[curindex][0]);
		}
		
		for (r=1;r<=R;r++)
		{
			if (n[r-1]!=0)
			{
			n[r-1]--; 
			int index_1r=popindex(n,R,planesizes); 
			n[r-1]++; // we have now the index of population N-1r

			mpf_t s; mpf_init(s); mpf_set_d(s,(double) qnm->L[m-1][r-1]);
			mpf_mul(s,s,G[index_1r][m_idx]);//qnm->L[m-1][r-1]);
			mpf_add(G[curindex][m_idx],G[curindex][m_idx],s);
			mpf_clear(s);
			}
		}
	}
	}
	while(!nextpop(n,N,R));
/*	for (r=1;r<=R;r++)
	{		
		if (N[r-1]>0)
		{
		N[r-1]--; 
		int index_1r=popindex(N,R,planesizes); 
		N[r-1]++; // we have now the index of population N-1r
		if (index_1r!=0)	
			{
				mpf_set(X[r-1],G[index_1r][M-1]);
				mpf_div(X[r-1],X[r-1],G[curindex][M-1]);
			}		
		else
			{
				mpf_set_d(X[r-1], (double) 1.0);
				mpf_div(X[r-1],X[r-1],G[curindex][M-1]);
			}
		}
		else
		{
			mpf_set_d(X[r-1],0.0);
		}
	}
*/
	// Return the normalizing constant G(N)
	if (g != NULL && g[0] != NULL) {
		int final_idx = hasZ ? M : M-1;
		mpf_set(g[0], G[curindex][final_idx]);
	}
	
	// Create Gk variables to store separate results for each k
	mpf_t** Gk=(mpf_t**) calloc(planesizes[0],sizeof(mpf_t*));
	for (t=1;t<=planesizes[0];t++)
	{
		Gk[t-1]=(mpf_t*) calloc(M,sizeof(mpf_t));
		for (m=1;m<=M;m++)
		{
			mpf_init(Gk[t-1][m-1]);
			mpf_set_d(Gk[t-1][m-1], (double) 0.0);
		}
	}
	for (m=1;m<=M;m++)
	{
		mpf_set_d(Gk[0][m-1], (double) 1.0);
	}
	
	// Reset population counter for independent convolutions
	int *n_k=(int*)initpop(R);
	int curindex_k=-1;
	
	// M independent convolutions using G values from previous loop
	for (int k=1; k<=M; k++)
	{
		// Reset for each k
		curindex_k=-1;
		for (t=1;t<=planesizes[0];t++)
		{
			for (m=1;m<=M;m++)
			{
				mpf_set_d(Gk[t-1][m-1], (double) 0.0);
			}
		}
		for (m=1;m<=M;m++)
		{
			mpf_set_d(Gk[0][m-1], (double) 1.0);
		}
		
		// Reset population for this k
		for (r=1;r<=R;r++)
		{
			n_k[r-1] = 0;
		}
		
		do
		{
			curindex_k++;
			
			int m_idx_k = hasZ ? k : k-1;
			
			if (k>1)
			{
				mpf_set(Gk[curindex_k][k-1],Gk[curindex_k][k-2]);
			}
			else if (hasZ)
			{
				mpf_set_d(Gk[curindex_k][k-1], (double) 1.0);
			}
			
			for (r=1;r<=R;r++)
			{
				if (n_k[r-1]!=0)
				{
					n_k[r-1]--; 
					int index_1r_k=popindex(n_k,R,planesizes); 
					n_k[r-1]++;
					
					mpf_t s_k; mpf_init(s_k); mpf_set_d(s_k,(double) qnm->L[k-1][r-1]);
					// Use G value from original computation at the end of loop at row 82
					mpf_mul(s_k,s_k,G[index_1r_k][m_idx_k]);
					mpf_add(Gk[curindex_k][k-1],Gk[curindex_k][k-1],s_k);
					mpf_clear(s_k);
				}
			}
		}
		while(!nextpop(n_k,N,R));
	}
	
	// Clean up Gk memory
	for (t=1;t<=planesizes[0];t++)
	{
		for (m=1;m<=M;m++)
		{
			mpf_clear(Gk[t-1][m-1]);
		}
		free(Gk[t-1]);
	}
	free(Gk);
	free(n_k);
	
	// Free memory
	for (t=1;t<=planesizes[0];t++)
	{
		for (m=1;m<=Mz;m++)
		{
			mpf_clear(G[t-1][m-1]);
		}
		free(G[t-1]);
	}
	free(G);
	free(n);
	free(planesizes);
	
	return;
}

