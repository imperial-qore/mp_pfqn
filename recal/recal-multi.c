#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gmpla.h"
#include "util.h"

void recal_multi(qnmodel* qn, double *Gex)
{
	int r;
	int M=qn->M;
	int R=qn->R;
	int Ntot=0;
	for (r=1;r<=R;r++)
		Ntot+=qn->N[r-1];
  	double* G, *G_1;
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

	G_1=(double*) calloc((long int)nck(Ntot+Mz-1,Ntot),sizeof(double));
	for (i=1;i<=nck(Ntot+Mz-1,Ntot);i++)
		G_1[i-1]=1;

	G=(double*) calloc((long int)nck(Ntot-1+Mz,Ntot),sizeof(double));

	int nr, n=0;
	int *m=calloc(Mz,sizeof(int));
	int *mZ=calloc(M,sizeof(int)); // For extracting first M elements when hasZ

	for (r=1;r<=R;r++) //for all classes
	{
		for (nr=1;nr<=qn->N[r-1];nr++)
		{
			n++;
			int** I_1=multichoose(Mz,(Ntot+1)-n);
			int** I=multichoose(Mz,(Ntot+1)-(n+1));
			for (i=1;i<=nck(Mz+(Ntot+1)-(n+1)-1,(Ntot+1)-(n+1));i++)
			{
				for (j=1; j<=Mz; j++)
					m[j-1]=I[i-1][j-1];		

				G[i-1]=0.0;
				
				if (hasZ && qn->Z[r-1] > 0) {
					// Add think time term: Z[r] * G_1[mzIndex] / nr
					// Extract first M elements for matching
					for (j=0; j<M; j++)
						mZ[j] = m[j];
					
					// Create temporary I_1 subset with only first M columns
					int** I_1_subset = calloc(nck(Mz+Ntot-n,Ntot-n+1), sizeof(int*));
					for (int k=0; k<nck(Mz+Ntot-n,Ntot-n+1); k++) {
						I_1_subset[k] = calloc(M, sizeof(int));
						for (int l=0; l<M; l++)
							I_1_subset[k][l] = I_1[k][l];
					}
					
					int mzIndex = int_matmatchrow(I_1_subset, nck(Mz+Ntot-n,Ntot-n+1), M, mZ);
					G[i-1] += ((double)qn->Z[r-1]) * G_1[mzIndex] / nr;
					
					// Free temporary subset
					for (int k=0; k<nck(Mz+Ntot-n,Ntot-n+1); k++)
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
					G[i-1]+=(m[j-1]+qn->mi[j-1]-1)*(qn->L[j-1][r-1])*(G_1[matchIndex])/nr;
					m[j-1]--;
				}
			}
			for (i=1;i<=nck(Mz+(Ntot+1)-(n+1)-1,(Ntot+1)-(n+1));i++)
				G_1[i-1]=G[i-1];
			free(I_1);
			I_1=I;
		}
	}
	*Gex=G[0];
	free(G);
	free(m);
	free(mZ);
	return;
}
