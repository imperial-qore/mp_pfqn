#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include "util.h"

double mva_multi(qnmodel* qn, mpq_t *X, mpq_t **Q, mpq_t G)
{
	int R=qn->R;
	int r,t,m;
	int M=qn->M;

	int *n=(int*)initpop(R);
	if (n == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory in mva_multi\n");
		return -1.0;
	}
  	int *planesizes=getplanesizes(qn->N,R); // planesize[k]=\prod_{j=R..k} (Nj+1)
	if (planesizes == NULL) {
		free(n);
		fprintf(stderr, "Error: Failed to allocate memory in mva_multi\n");
		return -1.0;
	}
	// init: the storage space using exact arithmetic
	int planesize0 = planesizes[0]; // Save value before potential free
	mpq_t** q=(mpq_t**) calloc(planesize0,sizeof(mpq_t*)); // the recursion is performed on N, so we keep the M queues of population N-1
	if (q == NULL) {
		free(n);
		free(planesizes);
		fprintf(stderr, "Error: Failed to allocate memory for q array (%d elements). Model may be too large.\n", planesize0);
		return -1.0;
	}

	for (t=1;t<=planesize0;t++) {
  		q[t-1]=(mpq_t*) calloc(M,sizeof(mpq_t));
		if (q[t-1] == NULL) {
			// Clean up previously allocated memory
			for (int i=0;i<t-1;i++) {
				for (m=0;m<M;m++)
					mpq_clear(q[i][m]);
				free(q[i]);
			}
			free(q);
			free(n);
			free(planesizes);
			fprintf(stderr, "Error: Failed to allocate memory for q[%d] (%d elements). Model may be too large.\n", t-1, M);
			return -1.0;
		}
		for (m=0;m<M;m++)
			mpq_init(q[t-1][m]);
	}
 	// end 
	int curindex=-1;
	
	// Initialize G = 1 at population N = [0,0,...,0] using exact arithmetic
	mpq_t G_accumulator;
	mpq_init(G_accumulator);
	mpq_set_ui(G_accumulator, 1, 1); // G starts at 1
	
	do	
	{
		curindex++; // current population index
		
		for (r=1;r<=R;r++)
		{
			if (n[r-1]!=0)
			{
			n[r-1]--; 
			int index_1r=popindex(n,R,planesizes); 
			n[r-1]++; // we have now the index of population N-1r
			
			// Compute Cr = Z[r] + sum(L[m][r] * (mi[m] + q[N-1r][m])) using exact arithmetic
			mpq_t Cr;
			mpq_init(Cr);
			mpq_set_z(Cr, qn->Z[r-1]);
			
			for (m=1;m<=M;m++)
			{
				mpq_t term;
				mpq_init(term);
				mpq_set_ui(term, qn->L[m-1][r-1] * qn->mi[m-1], 1);
				mpq_add(Cr, Cr, term);
				
				mpq_t L_times_q;
				mpq_init(L_times_q);
				mpq_set_ui(L_times_q, qn->L[m-1][r-1], 1);
				mpq_mul(L_times_q, L_times_q, q[index_1r][m-1]);
				mpq_add(Cr, Cr, L_times_q);
				
				mpq_clear(term);
				mpq_clear(L_times_q);
			}
			
			// X[r-1] = n[r-1] / Cr
			mpq_set_ui(X[r-1], n[r-1], 1);
			mpq_div(X[r-1], X[r-1], Cr);
			
			// Update queue lengths: q[curindex][m-1] += X[r-1] * L[m-1][r-1] * (mi[m-1] + q[index_1r][m-1])
			for (m=1;m<=M;m++) {
				mpq_t update_term;
				mpq_init(update_term);
				
				// L[m-1][r-1] * (mi[m-1] + q[index_1r][m-1])
				mpq_t mi_plus_q;
				mpq_init(mi_plus_q);
				mpq_set_ui(mi_plus_q, qn->mi[m-1], 1);
				mpq_add(mi_plus_q, mi_plus_q, q[index_1r][m-1]);
				
				mpq_set_ui(update_term, qn->L[m-1][r-1], 1);
				mpq_mul(update_term, update_term, mi_plus_q);
				mpq_mul(update_term, update_term, X[r-1]);
				
				mpq_add(q[curindex][m-1], q[curindex][m-1], update_term);
				
				mpq_clear(update_term);
				mpq_clear(mi_plus_q);
			}
			
			mpq_clear(Cr);
			}
		}
		
		// Find all non-zero indices and get the last valid one (matching Kotlin logic)
		int nonzero_indices[R];
		int nonzero_count = 0;
		
		// Find non-zero indices
		for (r = 0; r < R; r++) {
			if (n[r] > 0) {
				nonzero_indices[nonzero_count++] = r;
			}
		}
		
		if (nonzero_count > 0) {
			// Get the last non-zero index (matching Kotlin: last_nnz)
			int last_nnz = nonzero_indices[nonzero_count - 1];
			
			// Sum from 0 to last_nnz-1 (exclusive of last_nnz)
			int sumn = 0, sumN = 0;
			for (r = 0; r < last_nnz; r++) {
				sumn += n[r];
				sumN += qn->N[r];
			}
			
			// Sum from last_nnz+1 to R-1
			int sumnprime = 0;
			for (r = last_nnz + 1; r < R; r++) {
				sumnprime += n[r];
			}
			
			// Check class completion condition
			if (sumn == sumN && sumnprime == 0) {
				// Update G /= X[last_nnz] using exact arithmetic
				mpq_div(G_accumulator, G_accumulator, X[last_nnz]);
			}
		}
	}
	while(!nextpop(n,qn->N,R));
	
	// Set final G using exact arithmetic
	mpq_set(G, G_accumulator);
	mpq_clear(G_accumulator);
	
	// Compute final queue lengths for the full population using exact arithmetic
	for (m=1;m<=M;m++) {
		for (r=1;r<=R;r++) {
			// Get population N-1r (one less customer of class r)
			int *n_minus_1r = (int*)malloc(R * sizeof(int));
			for (int i=0; i<R; i++) {
				n_minus_1r[i] = qn->N[i];
			}
			n_minus_1r[r-1]--; // subtract one customer from class r
			
			int index_1r = popindex(n_minus_1r, R, planesizes);
			
			// Q[m-1][r-1] = X[r-1] * L[m-1][r-1] * (mi[m-1] + q[index_1r][m-1])
			mpq_t mi_plus_q;
			mpq_init(mi_plus_q);
			mpq_set_ui(mi_plus_q, qn->mi[m-1], 1);
			mpq_add(mi_plus_q, mi_plus_q, q[index_1r][m-1]);
			
			mpq_set_ui(Q[m-1][r-1], qn->L[m-1][r-1], 1);
			mpq_mul(Q[m-1][r-1], Q[m-1][r-1], mi_plus_q);
			mpq_mul(Q[m-1][r-1], Q[m-1][r-1], X[r-1]);
			
			mpq_clear(mi_plus_q);
			free(n_minus_1r);
		}
	}
	
	// Free memory
	for (t=1;t<=planesize0;t++) {
		for (m=0;m<M;m++)
			mpq_clear(q[t-1][m]);
		free(q[t-1]);
	}
	free(q);
	free(n);
	free(planesizes);
	
	// Compute logG from the final G value for return
	mpf_t G_mpf;
	mpf_init(G_mpf);
	mpf_set_q(G_mpf, G);
	double logG_final = log(mpf_get_d(G_mpf));
	mpf_clear(G_mpf);
	
	return logG_final;
}
