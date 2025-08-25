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
					mpz_mul(num, num, qnm->Z[r-1]);
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
			mpq_set_z(s, qnm->L[m-1][r-1]);
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
	
	// Compute Q[k][r] using checkpointed state and additional convolutions
	// After terminating the computation, save state of G variables
	// For each station k, resume from checkpoint and perform additional convolution
	// with extra queue type M+1 having identical demand as station k
	if (Q != NULL) {
		// Create checkpoint of current G state - deep copy entire G array
		mpq_t** G_checkpoint = (mpq_t**) calloc(planesize0, sizeof(mpq_t*));
		if (G_checkpoint == NULL) {
			fprintf(stderr, "Error: Failed to allocate memory for G checkpoint\n");
			return;
		}
		
		// Deep copy G array for checkpoint  
		for (t=1;t<=planesize0;t++) {
			G_checkpoint[t-1] = (mpq_t*) calloc(Mz+1, sizeof(mpq_t)); // +1 for extra station
			if (G_checkpoint[t-1] == NULL) {
				// Clean up partial checkpoint
				for (int j=0;j<t-1;j++) {
					for (m=0;m<Mz+1;m++) {
						mpq_clear(G_checkpoint[j][m]);
					}
					free(G_checkpoint[j]);
				}
				free(G_checkpoint);
				return;
			}
			for (m=0;m<Mz+1;m++) {
				mpq_init(G_checkpoint[t-1][m]);
				if (m < Mz) {
					mpq_set(G_checkpoint[t-1][m], G[t-1][m]);
				} else {
					mpq_set_ui(G_checkpoint[t-1][m], 0, 1); // Initialize extra station
				}
			}
		}
		
		// For each station k, compute Q[k][r] = mi[k] * G^k(N-er) / G(N)
		for (int k = 0; k < M; k++) {
			for (int r = 0; r < R; r++) {
				if (N[r] > 0) {
					// Create population N-er
					N[r]--;
					
					// Restore from checkpoint and extend G array temporarily  
					for (t=1;t<=planesize0;t++) {
						// Reallocate G[t-1] to have space for extra station
						mpq_t* temp_row = (mpq_t*) realloc(G[t-1], (Mz+1) * sizeof(mpq_t));
						if (temp_row == NULL) {
							N[r]++; // restore population
							continue;
						}
						G[t-1] = temp_row;
						
						// Initialize the new extra station slot
						mpq_init(G[t-1][Mz]);
						
						// Copy from checkpointed values
						for (m=0;m<Mz;m++) {
							mpq_set(G[t-1][m], G_checkpoint[t-1][m]);
						}
					}
					
					// Reset population iterator to beginning for the reduced population N-er
					int* n_iter = (int*)initpop(R);  
					if (n_iter == NULL) {
						N[r]++; // restore
						continue;
					}
					
					// Perform additional convolution for extra station with station k's demands
					int iter_index = -1;
					do {
						iter_index++;
						
						// Check if current population state n_iter is valid for N-er
						bool valid_state = true;
						for (int class_idx = 0; class_idx < R; class_idx++) {
							if (n_iter[class_idx] > N[class_idx]) {
								valid_state = false;
								break;
							}
						}
						if (!valid_state) continue;
						
						int current_pop_index = popindex(n_iter, R, planesizes);
						if (current_pop_index < 0 || current_pop_index >= planesize0) continue;
						
						// Initialize G for extra station from final original station
						mpq_set(G[current_pop_index][Mz], G[current_pop_index][Mz-1]);
						
						// Apply convolution step for extra station using station k's demands
						for (int class_idx = 1; class_idx <= R; class_idx++) {
							if (n_iter[class_idx-1] != 0) {
								n_iter[class_idx-1]--;
								int reduced_pop_index = popindex(n_iter, R, planesizes);
								n_iter[class_idx-1]++; // restore
								
								if (reduced_pop_index >= 0 && reduced_pop_index < planesize0) {
									mpq_t contribution;
									mpq_init(contribution);
									// Use station k's demand L[k][class_idx-1]
									mpq_set_z(contribution, qnm->L[k][class_idx-1]);
									mpq_mul(contribution, contribution, G[reduced_pop_index][Mz]);
									mpq_add(G[current_pop_index][Mz], G[current_pop_index][Mz], contribution);
									mpq_clear(contribution);
								}
							}
						}
						
					} while (!nextpop(n_iter, N, R));
					
					free(n_iter);
					
					// Get G^k(N-er) from the extra station at population N-er
					int final_pop_index = popindex(N, R, planesizes);  // N is still N-er here
					mpq_t G_k_N_minus_er;
					mpq_init(G_k_N_minus_er);
					
					if (final_pop_index >= 0 && final_pop_index < planesize0) {
						mpq_set(G_k_N_minus_er, G[final_pop_index][Mz]);
					} else {
						mpq_set_ui(G_k_N_minus_er, 0, 1);
					}
					
					// Restore population N
					N[r]++;
					
					// Compute Q[k][r] = L[k][r] * G^k(N-er) / G(N)  
					mpq_t q_exact;
					mpq_init(q_exact);
					
					// L[k][r] * G^k(N-er)  
					mpq_set_z(q_exact, qnm->L[k][r]);
					mpq_mul(q_exact, q_exact, G_k_N_minus_er);
					
					// Divide by G(N)
					if (g != NULL && g[0] != NULL && mpq_cmp_ui(g[0], 0, 1) != 0) {
						mpq_div(q_exact, q_exact, g[0]);
					} else {
						mpq_set_ui(q_exact, 0, 1);
					}
					
					// Convert to mpf_t
					mpf_set_q(Q[k][r], q_exact);
					
					mpq_clear(q_exact);
					mpq_clear(G_k_N_minus_er);
					
					// Clean up the extra station from G array
					for (t=1;t<=planesize0;t++) {
						mpq_clear(G[t-1][Mz]);
						// Resize back to original size
						mpq_t* temp_row = (mpq_t*) realloc(G[t-1], Mz * sizeof(mpq_t));
						if (temp_row != NULL) {
							G[t-1] = temp_row;
						}
					}
					
				} else {
					// If N[r] = 0, Q[k][r] = 0
					mpf_set_ui(Q[k][r], 0);
				}
			}
		}
		
		// Clean up checkpoint
		for (t=1;t<=planesize0;t++) {
			for (m=0;m<Mz+1;m++) {
				mpq_clear(G_checkpoint[t-1][m]);
			}
			free(G_checkpoint[t-1]);
		}
		free(G_checkpoint);
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
