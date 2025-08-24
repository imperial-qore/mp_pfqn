#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gmpla.h"
#include "util.h"

void recal_multi_exact(qnmodel* qn, mpq_t Gex)
{
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
		if (qn->Z[r-1] > 0) {
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
				return;
			}
			for (i=1;i<=nck(Mz+(Ntot+1)-(n+1)-1,(Ntot+1)-(n+1));i++)
			{
				for (j=1; j<=Mz; j++)
					m[j-1]=I[i-1][j-1];		

				mpq_set_ui(G[i-1], 0, 1);
				
				if (hasZ && qn->Z[r-1] > 0) {
					// Add think time term: Z[r] * G_1[mzIndex] / nr
					// Extract first M elements for matching
					for (j=0; j<M; j++)
						mZ[j] = m[j];
					
					// Create temporary I_1 subset with only first M columns
					int subset_size = nck(Mz+Ntot-n,Ntot-n+1);
					int** I_1_subset = calloc(subset_size, sizeof(int*));
					if (I_1_subset == NULL) {
						fprintf(stderr, "Error: Failed to allocate memory in recal. Model may be too large.\n");
						// Clean up and return
						for (i=1;i<=G_1_size;i++)
							mpq_clear(G_1[i-1]);
						free(G_1);
						for (i=1;i<=G_size;i++)
							mpq_clear(G[i-1]);
						free(G);
						free(m);
						free(mZ);
						free(I_1);
						free(I);
						return;
					}
					for (int k=0; k<subset_size; k++) {
						I_1_subset[k] = calloc(M, sizeof(int));
						if (I_1_subset[k] == NULL) {
							fprintf(stderr, "Error: Failed to allocate memory in recal. Model may be too large.\n");
							// Clean up
							for (int kk=0; kk<k; kk++)
								free(I_1_subset[kk]);
							free(I_1_subset);
							for (i=1;i<=G_1_size;i++)
								mpq_clear(G_1[i-1]);
							free(G_1);
							for (i=1;i<=G_size;i++)
								mpq_clear(G[i-1]);
							free(G);
							free(m);
							free(mZ);
							free(I_1);
							free(I);
							return;
						}
						for (int l=0; l<M; l++)
							I_1_subset[k][l] = I_1[k][l];
					}
					
					int mzIndex = int_matmatchrow(I_1_subset, nck(Mz+Ntot-n,Ntot-n+1), M, mZ);
					mpq_t term;
					mpq_init(term);
					mpq_set_ui(term, qn->Z[r-1], nr);
					mpq_mul(term, term, G_1[mzIndex]);
					mpq_add(G[i-1], G[i-1], term);
					mpq_clear(term);
					
					// Free temporary subset
					for (int k=0; k<subset_size; k++)
						free(I_1_subset[k]);
					free(I_1_subset);
				}
				
				for (j=1; j<=M; j++) // Only iterate over actual queues (not delay server)
				{
					m[j-1]++;
					int matchIndex;
					if (hasZ) {
						// Full match with M+1 elements
						matchIndex = int_matmatchrow(I_1, nck(Mz+Ntot-n,Ntot-n+1), Mz, m);
					} else {
						matchIndex = int_matmatchrow(I_1, nck(M+Ntot-n,Ntot-n+1), M, m);
					}
					
					mpq_t term;
					mpq_init(term);
					mpq_set_ui(term, (m[j-1]+qn->mi[j-1]-1)*(qn->L[j-1][r-1]), nr);
					mpq_mul(term, term, G_1[matchIndex]);
					mpq_add(G[i-1], G[i-1], term);
					mpq_clear(term);
					
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
	return;
}