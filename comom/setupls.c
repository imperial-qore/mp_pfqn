#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"
#include "profiling.h"

/* setuls - setup the linear system for the current class */
LS* setupls(combsrep *Dn, qnmodel* qnm, int *N, int r, double setup_start, int show_progress)
{
	int d,k,s;
	int M=qnm->M;
	int R=r;
	mpz_t** L=qnm->L;
	mpz_t* Z=qnm->Z;
	int* mi=qnm->mi;

	LS* ls= (LS*) calloc(1,sizeof(LS));

	long int cardG = nck(qnm->M+r-1,qnm->M);
	long int cardGk = cardG*qnm->M;

	long int cardCE = nck(qnm->M+r-2,qnm->M-1)*M;
	long int cardPC = nck(qnm->M+r-2,qnm->M-1)*(r-1);
 
 
	/* ginit */
	ls->A11=(mpq_mat_t)mpq_matzeros(cardGk,cardGk);
	ls->A12=(mpq_msp_t)mpq_msp(cardGk,cardG,cardCE+2*cardPC);
	ls->B1=(mpq_msp_t)mpq_msp(cardGk,cardGk+cardG,cardCE);
	ls->B2=(mpq_msp_t)mpq_msp(cardG,cardGk+cardG,cardG*(M+1));

	int * comb = calloc(sizeof(int),r);

	int row=0;
	struct rusage ruse;
	
	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum(Dn->combs[d-1],1-1,r-1)<M)
		{
			// Update progress on every iteration
			if (show_progress) {
				double current_time = CPUTIME;
				fprintf(stdout,"\r\033[K");
				int s;
				fprintf(stdout,"n=(");
				for (s=0;s<r-1;s++)
				    fprintf(stdout,"%d,",N[s]);
				fprintf(stdout,"%d) - %.6f s [Basis setup: %.6f s]",N[r-1], current_time - t0, current_time - setup_start);
				fflush(stdout);
			}
			for(s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			for (k=1; k<= M; k++)
			{
				/* add CE */
				row = row + 1;
				
				// Update progress on every inner loop iteration
				if (show_progress) {
					double current_time = CPUTIME;
					fprintf(stdout,"\r\033[K");
					int s;
					fprintf(stdout,"n=(");
					for (s=0;s<r-1;s++)
					    fprintf(stdout,"%d,",N[s]);
					fprintf(stdout,"%d) - %.6f s [Basis setup: %.6f s]",N[r-1], current_time - t0, current_time - setup_start);
					fflush(stdout);
				}
				mpq_matset_si(ls->A11,row-1,hash(Dn,comb,k)-1,1,1);
				mpq_mspset_si(ls->A12,row-1,hash(Dn,comb,0)-1-(int)nck(M+R-1,M)*M,-1,1);
				for (s=1;s<=r-1;s++)
				{
					comb[s-1] = comb[s-1] + 1;
					// Use mpz_t version to avoid overflow
					mpz_t neg_L;
					mpz_init(neg_L);
					mpz_neg(neg_L, L[k-1][s-1]);
					mpq_set_z(ls->A11[row-1][hash(Dn,comb,k)-1], neg_L);
					mpz_clear(neg_L);
					comb[s-1] = comb[s-1] - 1; 
					
				}
				mpq_mspset_z(ls->B1,row-1,hash(Dn,comb,k)-1,L[k-1][r-1],1);
			}
			for (s=1;s<=r-1;s++)
			{
				/* add PC */
				row = row + 1;
				mpq_mspset_si(ls->A12,row-1,hash(Dn,comb,0)-(int)nck(M+R-1,M)*M-1,N[s-1]-comb[s-1],1);	
				comb[s-1] = comb[s-1] + 1; /* oner(n,s) */
				// Use mpz_t version to avoid overflow
				mpz_t neg_Z;
				mpz_init(neg_Z);
				mpz_neg(neg_Z, Z[s-1]);
				mpq_mspset_z(ls->A12,row-1,hash(Dn,comb,0)-(int)nck(M+R-1,M)*M-1,neg_Z,1);
				mpz_clear(neg_Z);
				for (k=1; k<=M; k++) {
					// Calculate -mi[k-1]*L[k-1][s-1] using mpz_t
					mpz_t neg_mi_L;
					mpz_init(neg_mi_L);
					mpz_mul_si(neg_mi_L, L[k-1][s-1], -mi[k-1]);
					mpq_set_z(ls->A11[row-1][hash(Dn,comb,k)-1], neg_mi_L);
					mpz_clear(neg_mi_L);
				}
				comb[s-1] = comb[s-1] - 1; /* oner(n,s) */
			}
			
		}
	
	row = 0;
	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum((int*)Dn->combs[d-1],r-1,R-1-1)<=0)
		{
			// Update progress on every iteration in second main loop
			if (show_progress) {
				double current_time = CPUTIME;
				fprintf(stdout,"\r\033[K");
				int s;
				fprintf(stdout,"n=(");
				for (s=0;s<r-1;s++)
				    fprintf(stdout,"%d,",N[s]);
				fprintf(stdout,"%d) - %.6f s [Basis setup: %.6f s]",N[r-1], current_time - t0, current_time - setup_start);
				fflush(stdout);
			}
			row = row + 1;
			for(s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			mpq_mspset_z(ls->B2,row-1,hash(Dn,comb,0)-1,Z[r-1],1);	
			for (k=1;k<=M;k++) {
				// Calculate mi[k-1]*L[k-1][r-1] using mpz_t
				mpz_t mi_L;
				mpz_init(mi_L);
				mpz_mul_si(mi_L, L[k-1][r-1], mi[k-1]);
				mpq_mspset_z(ls->B2,row-1,hash(Dn,comb,k)-1,mi_L,1);
				mpz_clear(mi_L);
			}	
		}

	return ls;
}

