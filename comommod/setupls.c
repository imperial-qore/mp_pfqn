#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"

/* setupls - build the CoMoM linear system for class r in Z/p.
 *
 * A transliteration of comom/setupls.c: same rows in the same order, same
 * columns, same coefficients, with mpq_t replaced by a residue in a
 * machine word.  The progress reporting is dropped, since an image is not
 * interactive and may be running on one of several threads. */
LS* setupls(combsrep *Dn, qnmodel* qnm, int *N, int r)
{
	int d,k,s;
	int M=qnm->M;
	int R=r;
	mpz_t** L=qnm->L;
	mpz_t* Z=qnm->Z;
	int* mi=qnm->mi;
	int row=0;
	int* comb;

	LS* ls = (LS*) calloc(1,sizeof(LS));

	long int cardG  = nck(qnm->M+r-1,qnm->M);
	long int cardGk = cardG*qnm->M;
	long int cardCE = nck(qnm->M+r-2,qnm->M-1)*M;
	long int cardPC = nck(qnm->M+r-2,qnm->M-1)*(r-1);

	ls->m = M;
	ls->r = r;
	ls->cardG  = cardG;
	ls->cardGk = cardGk;

	ls->A11 = (zp_mat_t) zp_matzeros(cardGk,cardGk);
	ls->A12 = (zp_msp_t) zp_msp(cardGk,cardG,cardCE+2*cardPC);
	ls->B1  = (zp_msp_t) zp_msp(cardGk,cardGk+cardG,cardCE);
	ls->B2  = (zp_msp_t) zp_msp(cardG,cardGk+cardG,cardG*(M+1));

	comb = (int*) calloc(r,sizeof(int));

	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum(Dn->combs[d-1],1-1,r-1)<M)
		{
			for (s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			for (k=1;k<=M;k++)
			{
				/* add CE */
				row = row + 1;
				zp_matset_si(ls->A11,row-1,hash(Dn,comb,k)-1,1,1);
				zp_mspset_si(ls->A12,row-1,hash(Dn,comb,0)-1-(int)nck(M+R-1,M)*M,-1,1);
				for (s=1;s<=r-1;s++)
				{
					mpz_t neg_L;
					comb[s-1] = comb[s-1] + 1;
					mpz_init(neg_L);
					mpz_neg(neg_L, L[k-1][s-1]);
					ls->A11[row-1][hash(Dn,comb,k)-1] = zp_from_z(neg_L);
					mpz_clear(neg_L);
					comb[s-1] = comb[s-1] - 1;
				}
				zp_mspset_z(ls->B1,row-1,hash(Dn,comb,k)-1,L[k-1][r-1],1);
			}
			for (s=1;s<=r-1;s++)
			{
				mpz_t neg_Z;
				/* add PC */
				row = row + 1;
				zp_mspset_si(ls->A12,row-1,hash(Dn,comb,0)-(int)nck(M+R-1,M)*M-1,N[s-1]-comb[s-1],1);
				comb[s-1] = comb[s-1] + 1; /* oner(n,s) */
				mpz_init(neg_Z);
				mpz_neg(neg_Z, Z[s-1]);
				zp_mspset_z(ls->A12,row-1,hash(Dn,comb,0)-(int)nck(M+R-1,M)*M-1,neg_Z,1);
				mpz_clear(neg_Z);
				for (k=1;k<=M;k++)
				{
					mpz_t neg_mi_L;
					mpz_init(neg_mi_L);
					mpz_mul_si(neg_mi_L, L[k-1][s-1], -mi[k-1]);
					ls->A11[row-1][hash(Dn,comb,k)-1] = zp_from_z(neg_mi_L);
					mpz_clear(neg_mi_L);
				}
				comb[s-1] = comb[s-1] - 1; /* oner(n,s) */
			}
		}

	row = 0;
	for (d=1;d<=Dn->card;d++)
		if (int_vecsubsum((int*)Dn->combs[d-1],r-1,R-1-1)<=0)
		{
			row = row + 1;
			for (s=1;s<=r;s++)
				comb[s-1]=Dn->combs[d-1][s-1];
			zp_mspset_z(ls->B2,row-1,hash(Dn,comb,0)-1,Z[r-1],1);
			for (k=1;k<=M;k++)
			{
				mpz_t mi_L;
				mpz_init(mi_L);
				mpz_mul_si(mi_L, L[k-1][r-1], mi[k-1]);
				zp_mspset_z(ls->B2,row-1,hash(Dn,comb,k)-1,mi_L,1);
				mpz_clear(mi_L);
			}
		}

	free(comb);
	return ls;
}

void freels(LS* ls)
{
	if (!ls) return;
	if (ls->A11) zp_matfree(ls->A11, (int)ls->cardGk);
	zp_mspfree(ls->A12);
	zp_mspfree(ls->B1);
	zp_mspfree(ls->B2);
	if (ls->btf) btf_free(ls->btf);
	if (ls->A11_lu) zp_matfree(ls->A11_lu, (int)ls->cardGk);
	free(ls->A11_lu_indices);
	free(ls);
}
