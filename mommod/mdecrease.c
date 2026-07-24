#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include <zpla.h>
#include "mommod.h"

int mdecrease(qnmodel* qnm, zp_vec_t G, zp_vec_t Gk, zp_vec_t g, zp_vec_t gr, zp_t* res) /* G and Gk are passed by reference */
{
	zp_t tmp; zp_init(tmp);
	zp_t tmp2; zp_init(tmp2);
	int t,i,s,l,j;
	
	/* create a vector with all G and Gs, 1<=s<=qnm->R constants */
	Gk=(zp_vec_t)zp_vec(nck(qnm->M+qnm->R-1,qnm->R)*(qnm->R+1),0,1);
	for (t=0;t<nck(qnm->M+qnm->R-1,qnm->R);t++)
	{
		zp_set(Gk[t*(qnm->R+1)+qnm->R],gr[t*qnm->R]);
		for(s=0;s<qnm->R;s++)
			zp_set(Gk[t*(qnm->R+1)+s],g[t*qnm->R+s]);
	}
	int Ntot=0;
	for (s=0;s<qnm->R;s++)
		Ntot+=qnm->N[s];
	
	for (l=qnm->R;l>=1;l--)
	{
		G=(zp_vec_t)zp_vec(nck(qnm->M+l-2,l-1)*(qnm->R+1),0,1);		
		combsrep Ik;
		combsrep I;
		Ik.combs=(int**)sortbynnzpos((int**)multichoose(qnm->M,l),nck(qnm->M+l-1,l),qnm->M);
		I.combs=(int**)sortbynnzpos((int**)multichoose(qnm->M,l-1),nck(qnm->M+l-2,l-1),qnm->M);
		for (i=0;i<nck(qnm->M+l-2,l-1);i++)
		{
			/* compute G from G+k using the convolution expression */
			I.combs[i][0]++;
			t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
			I.combs[i][0]--;
			zp_set(G[i*(qnm->R+1)],Gk[t*(qnm->R+1)]);
			for (s=0;s<qnm->R;s++)
			{
				zp_set_z(tmp,qnm->L[0][s]);
				zp_mul(tmp,tmp,Gk[t*(qnm->R+1)+1+s]);		
				zp_sub(G[i*(qnm->R+1)],G[i*(qnm->R+1)],tmp);
			}

			/* compute Gs using the population contraints */
			for (s=1;s<=qnm->R;s++)
			{
				if (mpz_cmp_ui(qnm->Z[s-1], 0) == 0)
				{
					zp_set_si(G[i*(qnm->R+1)+s],0,1);
					int msum=0; // sum of multiplicities
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						zp_set_si(tmp,qnm->mi[j]+I.combs[i][j],1);
						zp_mul(tmp,tmp,Gk[t*(qnm->R+1)+s]);
						zp_add(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
						msum+=qnm->mi[j];
					}
					for (j=0;j<qnm->M;j++)
						msum+=I.combs[i][j];
					zp_set_si(tmp,Ntot+msum-1,1);
					zp_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);

				}
				else
				{
					zp_set_si(G[i*(qnm->R+1)+s],qnm->N[s-1],1);
					zp_set_z(tmp,qnm->Z[s-1]);
					zp_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
					zp_mul(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],G[i*(qnm->R+1)]);
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						zp_set_z(tmp,qnm->L[j][s-1]);
						zp_set_si(tmp2,qnm->mi[j]+I.combs[i][j],1);
						zp_mul(tmp2,tmp2,tmp);
						zp_set_z(tmp,qnm->Z[s-1]);
						zp_div(tmp2,tmp2,tmp);
						zp_mul(tmp2,tmp2,Gk[t*(qnm->R+1)+s]);
						zp_sub(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp2);
					}
				}
			}
		}
		if(l>1)	zp_vecdup(Gk,G,nck(qnm->M+l-2,l-1)*(qnm->R+1));
	}
	zp_clear(tmp);
	zp_clear(tmp2);
	exportres(qnm,G,Gk,res);
	return 0;
}

