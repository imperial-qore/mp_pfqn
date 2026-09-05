#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include <fpla.h>
#include "momf.h"

int mdecrease(qnmodel* qnm, fp_vec_t G, fp_vec_t Gk, fp_vec_t g, fp_vec_t gr, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor) /* G and Gk are passed by reference */
{
	fp_t tmp; fp_init(tmp);
	fp_t tmp2; fp_init(tmp2);
	int t,i,s,l,j;
	
	/* create a vector with all G and Gs, 1<=s<=qnm->R constants */
	Gk=(fp_vec_t)fp_vec(nck(qnm->M+qnm->R-1,qnm->R)*(qnm->R+1),0,1);
	for (t=0;t<nck(qnm->M+qnm->R-1,qnm->R);t++)
	{
		fp_set(Gk[t*(qnm->R+1)+qnm->R],gr[t*qnm->R]);
		for(s=0;s<qnm->R;s++)
			fp_set(Gk[t*(qnm->R+1)+s],g[t*qnm->R+s]);
	}
	int Ntot=0;
	for (s=0;s<qnm->R;s++)
		Ntot+=qnm->N[s];
	
	for (l=qnm->R;l>=1;l--)
	{
		G=(fp_vec_t)fp_vec(nck(qnm->M+l-2,l-1)*(qnm->R+1),0,1);		
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
			fp_set(G[i*(qnm->R+1)],Gk[t*(qnm->R+1)]);
			for (s=0;s<qnm->R;s++)
			{
				fp_set_z(tmp,qnm->L[0][s]);
				fp_mul(tmp,tmp,Gk[t*(qnm->R+1)+1+s]);		
				fp_sub(G[i*(qnm->R+1)],G[i*(qnm->R+1)],tmp);
			}

			/* compute Gs using the population contraints */
			for (s=1;s<=qnm->R;s++)
			{
				if (mpz_cmp_ui(qnm->Z[s-1], 0) == 0)
				{
					fp_set_si(G[i*(qnm->R+1)+s],0,1);
					int msum=0; // sum of multiplicities
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						fp_set_si(tmp,qnm->mi[j]+I.combs[i][j],1);
						fp_mul(tmp,tmp,Gk[t*(qnm->R+1)+s]);
						fp_add(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
						msum+=qnm->mi[j];
					}
					for (j=0;j<qnm->M;j++)
						msum+=I.combs[i][j];
					fp_set_si(tmp,Ntot+msum-1,1);
					fp_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);

				}
				else
				{
					fp_set_si(G[i*(qnm->R+1)+s],qnm->N[s-1],1);
					fp_set_z(tmp,qnm->Z[s-1]);
					fp_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
					fp_mul(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],G[i*(qnm->R+1)]);
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						fp_set_z(tmp,qnm->L[j][s-1]);
						fp_set_si(tmp2,qnm->mi[j]+I.combs[i][j],1);
						fp_mul(tmp2,tmp2,tmp);
						fp_set_z(tmp,qnm->Z[s-1]);
						fp_div(tmp2,tmp2,tmp);
						fp_mul(tmp2,tmp2,Gk[t*(qnm->R+1)+s]);
						fp_sub(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp2);
					}
				}
			}
		}
		if(l>1)	fp_vecdup(Gk,G,nck(qnm->M+l-2,l-1)*(qnm->R+1));
	}
	fp_clear(tmp);
	fp_clear(tmp2);
	perfindices(qnm,G,Gk,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
	return 0;
}

