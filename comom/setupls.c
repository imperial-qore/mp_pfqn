#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"

/* setuls - setup the linear system for the current class */
LS* setupls(combsrep *Dn, qnmodel* qnm, int *N, int r)
{
	int d,k,s;
	int M=qnm->M;
	int R=r;
	int** L=qnm->L;
	int* Z=qnm->Z;
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
	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum(Dn->combs[d-1],1-1,r-1)<M)
		{
			for(s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			for (k=1; k<= M; k++)
			{
				/* add CE */
				row = row + 1;
				mpq_matset_si(ls->A11,row-1,hash(Dn,comb,k+1)-1,1,1);
				mpq_mspset_si(ls->A12,row-1,hash(Dn,comb,0+1)-1-(int)nck(M+R-1,M)*M,-1,1);
				for (s=1;s<=r-1;s++)
				{
					comb[s-1] = comb[s-1] + 1;
					mpq_matset_si(ls->A11,row-1,hash(Dn,comb,k+1)-1,-L[k-1][s-1],1);
					comb[s-1] = comb[s-1] - 1; 
					
				}
				mpq_mspset_si(ls->B1,row-1,hash(Dn,comb,k+1)-1,L[k-1][r-1],1);
			}
			for (s=1;s<=r-1;s++)
			{
				/* add PC */
				row = row + 1;
				mpq_mspset_si(ls->A12,row-1,hash(Dn,comb,0+1)-(int)nck(M+R-1,M)*M-1,N[s-1]-comb[s-1],1);	
				comb[s-1] = comb[s-1] + 1; /* oner(n,s) */
				mpq_mspset_si(ls->A12,row-1,hash(Dn,comb,0+1)-(int)nck(M+R-1,M)*M-1,-Z[s-1],1);
				for (k=1; k<=M; k++)
					mpq_matset_si(ls->A11,row-1,hash(Dn,comb,k+1)-1,-mi[k-1]*L[k-1][s-1],1);
				comb[s-1] = comb[s-1] - 1; /* oner(n,s) */
			}
			
		}
	
	row = 0;
	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum((int*)Dn->combs[d-1],r-1,R-1-1)<=0)
		{
			row = row + 1;
			for(s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			mpq_mspset_si(ls->B2,row-1,hash(Dn,comb,0+1)-1,Z[r-1],1);	
			for (k=1;k<=M;k++)
				mpq_mspset_si(ls->B2,row-1,hash(Dn,comb,k+1)-1,mi[k-1]*L[k-1][r-1],1);	
		}

	return ls;
}

