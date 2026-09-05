#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "gmpla.h"
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

// Function to compute marginal normalizing constant G(N-e_r) for throughput calculation
void recal_marginal_exact(qnmodel* qn, int class_to_reduce, mpq_t G_marginal);

// Function to compute G^k(N-e_r) for queue length calculation
// Reduces population in class_to_reduce by 1 and increases multiplicity of station k by 1
void recal_queue_marginal_exact(qnmodel* qn, int class_to_reduce, int station_k, mpq_t G_k_marginal);

// Direct rank of a composition in the multichoose(ncols, ksum) enumeration
// order (v[.][0]=0..ksum, then recurse on the rest). Replaces the previous
// O(nrows) linear scan (int_matmatchrow); the outer i-loop called that per
// row, so the old code was O(nrows^2) per step and hung on think-time and
// large models. This rank is O(ncols + ksum).
static long int multichoose_rank(int ncols, int ksum, int* row)
{
	long int idx = 0;
	int pos = 0;
	while (ncols > 1) {
		int c = row[pos];
		for (int i = 0; i < c; i++)
			idx += nck((ncols-1)+(ksum-i)-1, ksum-i);
		ksum -= c;
		ncols -= 1;
		pos += 1;
	}
	return idx;
}

void recal_multi_exact_internal(qnmodel* qn, mpq_t Gex, bool show_progress)
{
	struct rusage ruse;
	double t_start = CPUTIME;
	int r;
	int M=qn->M;
	int R=qn->R;
	int Ntot=0;
	for (r=1;r<=R;r++)
		Ntot+=qn->N[r-1];
  	mpq_t* G, *G_1;
	int i, j;

	// Check if any think times are non-zero
	int hasZ = 0;
	for (r=1; r<=R; r++) {
		if (mpz_cmp_ui(qn->Z[r-1], 0) > 0) {
			hasZ = 1;
			break;
		}
	}

	// If think times are present, use M+1 stations (including delay server)
	int Mz = hasZ ? M + 1 : M;

	nckinit(1000000,1000000);

	long int G_1_size = nck(Ntot+Mz-1,Ntot);
	long int G_size = nck(Ntot-1+Mz,Ntot);
	
	if (G_1_size <= 0 || G_size <= 0) {
		fprintf(stderr, "Error: Invalid model dimensions for RECAL algorithm\n");
		return;
	}
	
	G_1=(mpq_t*) calloc((long int)G_1_size,sizeof(mpq_t));
	if (G_1 == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for G_1 array (%ld elements). Model may be too large.\n", G_1_size);
		return;
	}
	for (i=1;i<=G_1_size;i++) {
		mpq_init(G_1[i-1]);
		mpq_set_ui(G_1[i-1], 1, 1);
	}

	G=(mpq_t*) calloc((long int)G_size,sizeof(mpq_t));
	if (G == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for G array (%ld elements). Model may be too large.\n", G_size);
		for (i=1;i<=G_1_size;i++)
			mpq_clear(G_1[i-1]);
		free(G_1);
		return;
	}
	for (i=1;i<=G_size;i++) {
		mpq_init(G[i-1]);
	}

	int nr, n=0;
	int *pop_vector=calloc(R,sizeof(int)); // To track current population
	int *m=calloc(Mz,sizeof(int));
	if (m == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for m array (%d elements)\n", Mz);
		for (i=1;i<=G_1_size;i++)
			mpq_clear(G_1[i-1]);
		free(G_1);
		for (i=1;i<=G_size;i++)
			mpq_clear(G[i-1]);
		free(G);
		return;
	}
	int *mZ=calloc(M,sizeof(int)); // For extracting first M elements when hasZ
	if (mZ == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for mZ array (%d elements)\n", M);
		free(m);
		free(pop_vector);
		for (i=1;i<=G_1_size;i++)
			mpq_clear(G_1[i-1]);
		free(G_1);
		for (i=1;i<=G_size;i++)
			mpq_clear(G[i-1]);
		free(G);
		return;
	}

	for (r=1;r<=R;r++) //for all classes
	{
		for (nr=1;nr<=qn->N[r-1];nr++)
		{
			n++;
			pop_vector[r-1] = nr; // Update current population
			
			// Print progress: population vector and elapsed time
			if (show_progress) {
				// Calculate width needed for each population value
				int* widths = (int*) calloc(R, sizeof(int));
				for (int pr=0;pr<R;pr++) {
					int val = qn->N[pr];
					widths[pr] = 1;
					while (val >= 10) {
						widths[pr]++;
						val /= 10;
					}
				}
				
				fprintf(stderr, "\rn=(");
				for (int pr=0;pr<R;pr++) {
					fprintf(stderr, "%*d", widths[pr], pop_vector[pr]);
					if (pr<R-1) fprintf(stderr, ",");
				}
				double t_current = CPUTIME;
				fprintf(stderr, ") - Time: %.2f s  ", t_current - t_start);
				fflush(stderr);
				free(widths);
			}
			int** I_1=multichoose(Mz,(Ntot+1)-n);
			if (I_1 == NULL) {
				fprintf(stderr, "Error: Failed to allocate memory in multichoose. Model may be too large.\n");
				for (i=1;i<=G_1_size;i++)
					mpq_clear(G_1[i-1]);
				free(G_1);
				for (i=1;i<=G_size;i++)
					mpq_clear(G[i-1]);
				free(G);
				free(m);
				free(mZ);
				free(pop_vector);
				return;
			}
			int** I=multichoose(Mz,(Ntot+1)-(n+1));
			if (I == NULL) {
				fprintf(stderr, "Error: Failed to allocate memory in multichoose. Model may be too large.\n");
				free(I_1);
				for (i=1;i<=G_1_size;i++)
					mpq_clear(G_1[i-1]);
				free(G_1);
				for (i=1;i<=G_size;i++)
					mpq_clear(G[i-1]);
				free(G);
				free(m);
				free(mZ);
				free(pop_vector);
				return;
			}
			for (i=1;i<=nck(Mz+(Ntot+1)-(n+1)-1,(Ntot+1)-(n+1));i++)
			{
				for (j=1; j<=Mz; j++)
					m[j-1]=I[i-1][j-1];		

				mpq_set_ui(G[i-1], 0, 1);
				
				if (hasZ && mpz_cmp_ui(qn->Z[r-1], 0) > 0) {
					// Think time term: Z[r] * G_1[mzIndex] / nr, where mzIndex
					// is the configuration with one extra job at the delay
					// server (last station, column index M). The matching row
					// is unique, so rank it directly instead of scanning.
					m[M]++;
					long int mzIndex = multichoose_rank(Mz, (Ntot+1)-n, m);
					m[M]--;
					mpq_t term, nr_q;
					mpq_init(term);
					mpq_init(nr_q);
					mpq_set_z(term, qn->Z[r-1]);
					mpq_set_ui(nr_q, nr, 1);
					mpq_div(term, term, nr_q);
					mpq_mul(term, term, G_1[mzIndex]);
					mpq_add(G[i-1], G[i-1], term);
					mpq_clear(term);
					mpq_clear(nr_q);
				}
				
				for (j=1; j<=M; j++) // Only iterate over actual queues (not delay server)
				{
					m[j-1]++;
					long int matchIndex;
					if (hasZ) {
						matchIndex = multichoose_rank(Mz, (Ntot+1)-n, m);
					} else {
						matchIndex = multichoose_rank(M, (Ntot+1)-n, m);
					}
					
					mpq_t term, nr_q, L_term;
					mpq_init(term);
					mpq_init(nr_q);
					mpq_init(L_term);
					mpq_set_z(L_term, qn->L[j-1][r-1]);
					mpq_set_ui(term, m[j-1]+qn->mi[j-1]-1, 1);
					mpq_mul(term, term, L_term);
					mpq_set_ui(nr_q, nr, 1);
					mpq_div(term, term, nr_q);
					mpq_mul(term, term, G_1[matchIndex]);
					mpq_add(G[i-1], G[i-1], term);
					mpq_clear(term);
					mpq_clear(nr_q);
					mpq_clear(L_term);
					
					m[j-1]--;
				}
			}
			// Copy G to G_1
			for (i=1;i<=nck(Mz+(Ntot+1)-(n+1)-1,(Ntot+1)-(n+1));i++) {
				mpq_set(G_1[i-1], G[i-1]);
			}
			free(I_1);
			I_1=I;
		}
	}
	
	// Clear the progress line before performance metrics display
	if (show_progress) {
		fprintf(stderr, "\r%*s\r", 80, "");
	}
	
	mpq_set(Gex, G[0]);
	
	// Free memory
	for (i=1;i<=nck(Ntot+Mz-1,Ntot);i++)
		mpq_clear(G_1[i-1]);
	free(G_1);
	for (i=1;i<=nck(Ntot-1+Mz,Ntot);i++)
		mpq_clear(G[i-1]);
	free(G);
	free(m);
	free(mZ);
	free(pop_vector);
	
	return;
}

void recal_multi_exact(qnmodel* qn, mpq_t Gex)
{
	recal_multi_exact_internal(qn, Gex, true);
}

void recal_marginal_exact(qnmodel* qn, int class_to_reduce, mpq_t G_marginal)
{
	// Create a temporary model with reduced population for the specified class
	qnmodel temp_qn = *qn;
	
	// Allocate new arrays for the temporary model
	temp_qn.N = (int*)malloc(qn->R * sizeof(int));
	for (int i = 0; i < qn->R; i++) {
		temp_qn.N[i] = qn->N[i];
	}
	
	// Reduce population by 1 for the specified class
	if (temp_qn.N[class_to_reduce] > 0) {
		temp_qn.N[class_to_reduce]--;
	} else {
		// If population is already 0, return 0
		mpq_set_ui(G_marginal, 0, 1);
		free(temp_qn.N);
		return;
	}
	
	// Compute G(N-e_r) using the modified population
	// Show progress for marginal computation
	fprintf(stderr, "\rComputing X[%d] marginal...                    ", class_to_reduce + 1);
	fflush(stderr);
	recal_multi_exact_internal(&temp_qn, G_marginal, false);
	
	// Free temporary arrays
	free(temp_qn.N);
}

void recal_queue_marginal_exact(qnmodel* qn, int class_to_reduce, int station_k, mpq_t G_k_marginal)
{
	// Create a temporary model with reduced population and increased multiplicity
	qnmodel temp_qn = *qn;
	
	// Allocate new arrays for the temporary model
	temp_qn.N = (int*)malloc(qn->R * sizeof(int));
	temp_qn.mi = (int*)malloc(qn->M * sizeof(int));
	
	// Copy original arrays
	for (int i = 0; i < qn->R; i++) {
		temp_qn.N[i] = qn->N[i];
	}
	for (int i = 0; i < qn->M; i++) {
		temp_qn.mi[i] = qn->mi[i];
	}
	
	// Reduce population by 1 for the specified class
	if (temp_qn.N[class_to_reduce] > 0) {
		temp_qn.N[class_to_reduce]--;
	} else {
		// If population is already 0, return 0
		mpq_set_ui(G_k_marginal, 0, 1);
		free(temp_qn.N);
		free(temp_qn.mi);
		return;
	}
	
	// Increase multiplicity by 1 for the specified station
	temp_qn.mi[station_k]++;
	
	// Compute G^k(N-e_r) using the modified population and multiplicity
	// Show progress for queue marginal computation
	fprintf(stderr, "\rComputing Q[%d,%d] marginal...                    ", station_k + 1, class_to_reduce + 1);
	fflush(stderr);
	recal_multi_exact_internal(&temp_qn, G_k_marginal, false);
	
	// Free temporary arrays
	free(temp_qn.N);
	free(temp_qn.mi);
}