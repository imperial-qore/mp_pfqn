#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"

int mdecrease(qnmodel* qnm, mpq_vec_t G, mpq_vec_t Gk, mpq_vec_t g, mpq_vec_t gr, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, bool debug_full_precision, mpz_t scale_factor) /* G and Gk are passed by reference */
{
	mpq_t tmp; mpq_init(tmp);
	mpq_t tmp2; mpq_init(tmp2);
	int t,i,s,l,j;
	
	/* create a vector with all G and Gs, 1<=s<=qnm->R constants */
	Gk=(mpq_vec_t)mpq_vec(nck(qnm->M+qnm->R-1,qnm->R)*(qnm->R+1),0,1);
	for (t=0;t<nck(qnm->M+qnm->R-1,qnm->R);t++)
	{
		mpq_set(Gk[t*(qnm->R+1)+qnm->R],gr[t*qnm->R]);
		for(s=0;s<qnm->R;s++)
			mpq_set(Gk[t*(qnm->R+1)+s],g[t*qnm->R+s]);
	}
	int Ntot=0;
	for (s=0;s<qnm->R;s++)
		Ntot+=qnm->N[s];
	
	for (l=qnm->R;l>=1;l--)
	{
		G=(mpq_vec_t)mpq_vec(nck(qnm->M+l-2,l-1)*(qnm->R+1),0,1);		
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
			mpq_set(G[i*(qnm->R+1)],Gk[t*(qnm->R+1)]);
			for (s=0;s<qnm->R;s++)
			{
				mpq_set_z(tmp,qnm->L[0][s]);
				mpq_mul(tmp,tmp,Gk[t*(qnm->R+1)+1+s]);		
				mpq_sub(G[i*(qnm->R+1)],G[i*(qnm->R+1)],tmp);
			}

			/* compute Gs using the population contraints */
			for (s=1;s<=qnm->R;s++)
			{
				if (mpz_cmp_ui(qnm->Z[s-1], 0) == 0)
				{
					mpq_set_si(G[i*(qnm->R+1)+s],0,1);
					int msum=0; // sum of multiplicities
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						mpq_set_si(tmp,qnm->mi[j]+I.combs[i][j],1);
						mpq_mul(tmp,tmp,Gk[t*(qnm->R+1)+s]);
						mpq_add(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
						msum+=qnm->mi[j];
					}
					for (j=0;j<qnm->M;j++)
						msum+=I.combs[i][j];
					mpq_set_si(tmp,Ntot+msum-1,1);
					mpq_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);

				}
				else
				{
					mpq_set_si(G[i*(qnm->R+1)+s],qnm->N[s-1],1);
					mpq_set_z(tmp,qnm->Z[s-1]);
					mpq_div(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp);
					mpq_mul(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],G[i*(qnm->R+1)]);
					for (j=0;j<qnm->M;j++)
					{
						I.combs[i][j]++;
						t=int_matmatchrow(Ik.combs,nck(qnm->M+l-1,l),qnm->M,I.combs[i]);
						I.combs[i][j]--;
						mpq_set_z(tmp,qnm->L[j][s-1]);
						mpq_set_si(tmp2,qnm->mi[j]+I.combs[i][j],1);
						mpq_mul(tmp2,tmp2,tmp);
						mpq_set_z(tmp,qnm->Z[s-1]);
						mpq_div(tmp2,tmp2,tmp);
						mpq_mul(tmp2,tmp2,Gk[t*(qnm->R+1)+s]);
						mpq_sub(G[i*(qnm->R+1)+s],G[i*(qnm->R+1)+s],tmp2);
					}
				}
			}
		}
		if(l>1)	mpq_vecdup(Gk,G,nck(qnm->M+l-2,l-1)*(qnm->R+1));
	}
	mpq_clear(tmp);
	mpq_clear(tmp2);
	
	// Print all normalizing constants if -d option is used
	if (debug_output && !log_output && !throughput_output && !queue_output) {
		printf("\n========== Normalizing Constants (Debug) ==========\n");
		
		// Print last basis (current g vector)
		printf("Last basis (final population):\n");
		long int final_cardG = nck(qnm->M+qnm->R-1,qnm->R-1);
		long int final_cardGk = final_cardG*qnm->R;
		
		for (int i = 0; i < final_cardGk + final_cardG; i++) {
			if (debug_full_precision) {
				gmp_printf("g[%d] = %Qd\n", i, g[i]);
			} else {
				mpf_t fval_debug;
				mpf_init(fval_debug);
				mpf_set_q(fval_debug, g[i]);
				printf("g[%d] = %.15e\n", i, mpf_get_d(fval_debug));
				mpf_clear(fval_debug);
			}
		}
		
		// Print penultimate basis (gr vector)
		if (qnm->N[qnm->R-1] > 1) {
			printf("\nPenultimate basis (population N[R]-1):\n");
		} else {
			printf("\nPenultimate basis (population N[R-1], since N[R]=1):\n");
		}
		for (int i = 0; i < final_cardGk + final_cardG; i++) {
			if (debug_full_precision) {
				gmp_printf("gr[%d] = %Qd\n", i, gr[i]);
			} else {
				mpf_t fval_debug;
				mpf_init(fval_debug);
				mpf_set_q(fval_debug, gr[i]);
				printf("gr[%d] = %.15e\n", i, mpf_get_d(fval_debug));
				mpf_clear(fval_debug);
			}
		}
		
		printf("===================================================\n");
	}
	
	perfindices(qnm,G,Gk,verbose_output,log_output,normconst_output,normconst_g_output,throughput_output,queue_output,debug_output,bounds_output,scale_factor);
	return 0;
}

