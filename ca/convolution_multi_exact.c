#include <stdio.h>
#include <stdlib.h>
#include <math.h>		
#include <gmp.h>
#include "mom.h"
#include "popcycle.h"

void convolution_multi_exact(mpq_t *g, mpf_t *X, mpf_t **Q)
{
	int M=qnm->M;
	int R=qnm->R;
	int r,t,m;
	int*N=qnm->N;
	int *n=(int*)initpop(R);
	if (n == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory in convolution_multi_exact\n");
		return;
	}
  	int *planesizes=(int*) getplanesizes(N,R); // planesize[k]=\prod_{j=R..k} (Nj+1)
	if (planesizes == NULL) {
		free(n);
		fprintf(stderr, "Error: Failed to allocate memory in convolution_multi_exact\n");
		return;
	}
	
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
	
	// init: the storage space using exact rationals
	int planesize0 = planesizes[0]; // Save value before potential free
	mpq_t** G=(mpq_t**) calloc(planesize0,sizeof(mpq_t*));
	if (G == NULL) {
		free(n);
		free(planesizes);
		fprintf(stderr, "Error: Failed to allocate memory for G array (%d elements). Model may be too large.\n", planesize0);
		return;
	}
	
	for (t=1;t<=planesize0;t++)
	{
  		G[t-1]=(mpq_t*) calloc(Mz,sizeof(mpq_t));
		if (G[t-1] == NULL) {
			// Clean up previously allocated memory
			for (int i=0;i<t-1;i++) {
				for (m=1;m<=Mz;m++) {
					mpq_clear(G[i][m-1]);
				}
				free(G[i]);
			}
			free(G);
			free(n);
			free(planesizes);
			fprintf(stderr, "Error: Failed to allocate memory for G[%d] (%d elements). Model may be too large.\n", t-1, Mz);
			return;
		}
		for (m=1;m<=Mz;m++)
		{
  			mpq_init(G[t-1][m-1]);
			mpq_set_ui(G[t-1][m-1], 0, 1);
		}
	}
	for (m=1;m<=Mz;m++)
	{
		mpq_set_ui(G[0][m-1], 1, 1);
	}
 	// end 
	int curindex=-1;
	do	
	{
	curindex++; // current population index
	
	// If think times present, first compute G[0,n] based on Z (delay server)
	if (hasZ) {
		mpq_set_ui(G[curindex][0], 1, 1);
		for (r=1;r<=R;r++)
		{
			if (n[r-1]!=0 && qnm->Z[r-1] > 0)
			{
				// G[0,n] = prod(Z[r]^n[r] / n[r]!)
				mpq_t zr_pow; mpq_init(zr_pow);
				mpq_set_ui(zr_pow, 1, 1);
				
				// Compute Z[r]^n[r] / n[r]!
				mpz_t num, den;
				mpz_init(num);
				mpz_init(den);
				mpq_get_num(num, zr_pow);
				mpq_get_den(den, zr_pow);
				
				for (int i = 1; i <= n[r-1]; i++) {
					mpz_mul_ui(num, num, qnm->Z[r-1]);
					mpz_mul_ui(den, den, i);
				}
				
				mpq_set_num(zr_pow, num);
				mpq_set_den(zr_pow, den);
				mpq_canonicalize(zr_pow);
				
				mpz_clear(num);
				mpz_clear(den);
				
				mpq_mul(G[curindex][0], G[curindex][0], zr_pow);
				mpq_clear(zr_pow);
			}
		}
	}
	
	for (m=1;m<=M;m++)
	{
		int m_idx = hasZ ? m : m-1; // Adjust index when delay server is present
		
		if (m>1) // add norm const of the complementary network
		{
			mpq_set(G[curindex][m_idx],G[curindex][m_idx-1]);
		}
		else if (hasZ) // m==1 with think times, copy from delay server
		{
			mpq_set(G[curindex][m_idx],G[curindex][0]);
		}
		
		for (r=1;r<=R;r++)
		{
			if (n[r-1]!=0)
			{
			n[r-1]--; 
			int index_1r=popindex(n,R,planesizes); 
			n[r-1]++; // we have now the index of population N-1r

			mpq_t s; mpq_init(s);
			mpq_set_ui(s, qnm->L[m-1][r-1], 1);
			mpq_mul(s,s,G[index_1r][m_idx]);
			mpq_add(G[curindex][m_idx],G[curindex][m_idx],s);
			mpq_clear(s);
			}
		}
	}
	}
	while(!nextpop(n,N,R));
	
	// Return the normalizing constant G(N)
	if (g != NULL && g[0] != NULL) {
		int final_idx = hasZ ? M : M-1;
		mpq_set(g[0], G[curindex][final_idx]);
	}
	
	// Compute throughputs X[s] = G(N-1s)/G(N)
	if (X != NULL) {
		for (r=1; r<=R; r++) {
			// Create N-1r population vector
			N[r-1]--; 
			int index_1r = popindex(N, R, planesizes);
			N[r-1]++; // restore N
			
			if (index_1r >= 0 && g != NULL && g[0] != NULL) {
				// X[r-1] = G(N-1r) / G(N)
				mpq_t g_1r;
				mpq_init(g_1r);
				int final_idx = hasZ ? M : M-1;
				mpq_set(g_1r, G[index_1r][final_idx]);
				
				// Convert to mpf_t and compute ratio
				mpf_t g_1r_f, g_n_f;
				mpf_init(g_1r_f);
				mpf_init(g_n_f);
				mpf_set_q(g_1r_f, g_1r);
				mpf_set_q(g_n_f, g[0]);
				
				mpf_div(X[r-1], g_1r_f, g_n_f);
				
				mpq_clear(g_1r);
				mpf_clear(g_1r_f);
				mpf_clear(g_n_f);
			}
		}
	}
	
	// Compute Q[i][r] = L[i][r] * G^+i_r / G
	// where G^+i_r is G at population N with queue i multiplicity increased by 1
	if (Q != NULL) {
		
		for (int i = 0; i < M; i++) {
			// Compute G^+i at population N with queue i multiplicity increased by 1
			// Start from G at station i-1 and redo convolution for station i with mi+1
			
			mpq_t G_plus_i;
			mpq_init(G_plus_i);
			
			// Start with G from previous station (or delay server if hasZ)
			if (i == 0) {
				// First station: start from delay server or initial condition
				if (hasZ) {
					mpq_set(G_plus_i, G[curindex][0]); // From delay server
				} else {
					mpq_set_ui(G_plus_i, 1, 1); // Initial condition
				}
			} else {
				// Copy from previous station
				mpq_set(G_plus_i, G[curindex][hasZ ? i : i-1]);
			}
			
			// Now apply convolution for station i with multiplicity increased by 1
			// This means we apply the convolution step twice (original + 1 additional)
			for (r = 1; r <= R; r++) {
				if (N[r-1] != 0) {
					N[r-1]--; 
					int index_1r = popindex(N, R, planesizes); 
					N[r-1]++; // restore N
					
					if (index_1r >= 0) {
						mpq_t s; 
						mpq_init(s);
						mpq_set_ui(s, qnm->L[i][r-1], 1);
						
						// Get G^+i at the reduced population from previous station
						mpq_t G_prev;
						mpq_init(G_prev);
						if (i == 0) {
							if (hasZ) {
								mpq_set(G_prev, G[index_1r][0]); // From delay server
							} else {
								mpq_set_ui(G_prev, 1, 1); // Initial condition
							}
						} else {
							mpq_set(G_prev, G[index_1r][hasZ ? i : i-1]);
						}
						
						mpq_mul(s, s, G_prev);
						mpq_add(G_plus_i, G_plus_i, s);
						
						mpq_clear(s);
						mpq_clear(G_prev);
					}
				}
			}
			
			
			// Now compute Q[i][r] for each class r using this G^+i
			for (r = 0; r < R; r++) {
				if (g != NULL && g[0] != NULL) {
					mpq_t q_exact;
					mpq_init(q_exact);
					
					// The issue is that we're using G^+i at population N
					// But we need G^+i at population N-1r
					// Let's compute the correct value
					
					if (N[r] > 0) {
						// We need to compute G^+i(N-1r)
						mpq_t G_plus_i_N_minus_1;
						mpq_init(G_plus_i_N_minus_1);
						
						// Get index for population N-1r
						N[r]--;
						int index_1r = popindex(N, R, planesizes);
						
						if (index_1r >= 0) {
							// Start from G(N-1r, stations 0..i-1)
							if (i == 0) {
								if (hasZ) {
									mpq_set(G_plus_i_N_minus_1, G[index_1r][0]);
								} else {
									mpq_set_ui(G_plus_i_N_minus_1, 1, 1);
								}
							} else {
								int prev_idx = hasZ ? i : i-1;
								mpq_set(G_plus_i_N_minus_1, G[index_1r][prev_idx]);
							}
							
							// Apply convolution for station i with the extra application
							for (int s = 0; s < R; s++) {
								if (N[s] > 0) { // N is already reduced by 1 for class r
									N[s]--;
									int index_2 = popindex(N, R, planesizes);
									N[s]++;
									
									if (index_2 >= 0) {
										mpq_t contrib;
										mpq_init(contrib);
										mpq_set_ui(contrib, qnm->L[i][s], 1);
										
										if (i == 0) {
											// No previous station
										} else {
											int prev_idx = hasZ ? i : i-1;
											mpq_mul(contrib, contrib, G[index_2][prev_idx]);
										}
										
										mpq_add(G_plus_i_N_minus_1, G_plus_i_N_minus_1, contrib);
										mpq_clear(contrib);
									}
								}
							}
							
						}
						
						N[r]++; // Restore N
						
						// Q[i][r] = L[i][r] * G^+i(N-1r) / G(N)
						// But based on ground truth, for station 2 in single class, we need a factor of (N-1)
						mpq_set_ui(q_exact, qnm->L[i][r], 1);
						
						mpq_mul(q_exact, q_exact, G_plus_i_N_minus_1);
						
						mpq_div(q_exact, q_exact, g[0]);
						
						mpq_clear(G_plus_i_N_minus_1);
					} else {
						// If N[r] = 0, Q[i][r] = 0
						mpq_set_ui(q_exact, 0, 1);
					}
					
					// Convert to mpf_t
					mpf_set_q(Q[i][r], q_exact);
					
					mpq_clear(q_exact);
				}
			}
			
			mpq_clear(G_plus_i);
		}
	}
	
	// Free memory
	for (t=1;t<=planesize0;t++)
	{
		for (m=1;m<=Mz;m++)
		{
			mpq_clear(G[t-1][m-1]);
		}
		free(G[t-1]);
	}
	free(G);
	free(n);
	free(planesizes);
	
	return;
}
