#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"

/* setuls - setup the linear system for the current class */
LS* setupls(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r)
{
	return setupls_progress(L, n, Z, mi, m, r, 0.0, 0);
}

LS* setupls_progress(mpz_t** L, int* n, mpz_t *Z, int* mi, int m, int r, double setup_start, int show_progress)
{
	int i,j,k,s,h,t,w;
	combsrep Ik;
	combsrep I;
	LS* ls= (LS*) calloc(1,sizeof(LS));
	int *pcpos;
	int lastA12,lastD,lastB2r,lastpcpos;
	ls->m=m;
	ls->r=r;
	ls->H=MIN(m,r);

	int* maxL=(int*)int_vec(r,1);
//	for (s=0;s<r;s++)
//		for (i=0;i<m;i++)
//		maxL[s]=MIN(maxL[s],L[i][s]);
	
	/* compute Ik */
	Ik.n=m;
	Ik.k=r;
	Ik.card=nck(m+r-1,r);
	Ik.combs=(int**)sortbynnzpos((int**)multichoose(Ik.n,Ik.k),Ik.card,Ik.n);

	/* compute I */
	I.n=m;
	I.k=r-1;
	I.card=nck(m+r-2,r-1);
	I.combs=(int**)sortbynnzpos((int**)multichoose(I.n,I.k),I.card,I.n);
                
	/* initialize C matrix */
	for (h=1;h<=ls->H;h++) ls->numdiagblocks += nck(m,h); /* compute number of diagonal blocks */
	ls->C=(mpq_diagblock*) calloc(ls->numdiagblocks+1,sizeof(mpq_diagblock));
	int* last=(int*) int_vec(ls->numdiagblocks,0);
	int* d2diagblock=(int*) int_vec(nck(m+r-1,r),0);
	int* d2relpos=(int*) int_vec(nck(m+r-1,r),0);
	t=0; 
	for (h=1;h<=ls->H;h++) /* for all diagonal blocks */
	{
		for (k=0;k<nck(m,h);k++)
		{
		ls->C[t+k].diag=(mpq_t**) mpq_matzeros(nck(r-1,r-h)*r,nck(r-1,r-h)*r);
		ls->C[t+k].nondiag=(mpq_msp_t) mpq_msp(nck(r-1,r-h)*r,Ik.card*r,nck(r-1,r-h)*(r-h)*(m-h));
		}
		t+=nck(m,h);
	}
	/* initialize A12, D and B2r matrices */
	ls->A12=(mpq_msp_t)mpq_msp(nck(m+r-1,r)*r,nck(m+r-2,r-1)*r,nck(m+r-2,r-1)*m+2*(nck(m+r-1,r)*r-nck(m+r-2,r-1)*m));
	lastA12=0;
	ls->B1r=(mpq_msp_t)mpq_msp(nck(m+r-1,r)*r,nck(m+r-1,r)*r+nck(m+r-2,r-1)*r,nck(m+r-2,r-1)*m);
	lastD=0;
	ls->B2r=(mpq_msp_t)mpq_msp(nck(m+r-2,r-1)*r,nck(m+r-1,r)*r+nck(m+r-2,r-1)*r,nck(m+r-2,r-1)*r*(m+1));
	lastB2r=0;
	pcpos=(int*)int_vec(nck(m+r-2,r-1)*(r-1),0);
	lastpcpos=0;

	/* linear system setup */
	t=0;w=0;
	for (h=1;h<=ls->H;h++)
	{
	for (k=0;k<nck(m,h);k++)
	{
		last[t]=0;
		for (i=0;i<nck(r-1,r-h);i++)
		{
			int* d = Ik.combs[w];
			d2diagblock[w]=t;
			d2relpos[w]=i;
			for (j=0;j<m;j++)
			{
			 if (d[j]>0)
			 {
			 	mpq_set_si(ls->C[t].diag[last[t]][i*r+0],1,1);
				for (s=1;s<r;s++)
			 		{
						mpz_t neg_L;
						mpz_init(neg_L);
						mpz_neg(neg_L, L[j][s-1]);
						mpq_t divisor;
						mpq_init(divisor);
						mpq_set_ui(divisor, maxL[s-1], 1);
						mpq_set_z(ls->C[t].diag[last[t]][i*r+s], neg_L);
						mpq_div(ls->C[t].diag[last[t]][i*r+s], ls->C[t].diag[last[t]][i*r+s], divisor);
						mpz_clear(neg_L);
						mpq_clear(divisor);
					}
				d[j]--;
				mpq_mspset_si(ls->A12,lastA12,int_matmatchrow(I.combs,I.card,I.n,d)*r,-1,1);
				d[j]++;
				mpq_mspset_z(ls->B1r,lastD,w*r,L[j][r-1],maxL[r-1]);
				lastA12++;
				lastD++;
				last[t]++;
			 }
			}
			w++;
		} /* for i */
		for (j=0;j<nck(r-1,r-h)*(r-h);j++)
		{
			lastD++;
			pcpos[lastpcpos++]=lastA12++;
		}
		t++;
	} /* for k */
	} /* for h */
	w=0;
	lastpcpos=0;	
	for (i=0;i<I.card;i++)
	{
		int* d = I.combs[i];
		t=0;
		for (s=0;s<=r-1;s++)
			mpq_mspset_z(ls->B2r,lastB2r++,(Ik.card+i)*r+s,Z[r-1],maxL[r-1]);
		
		for (j=0;j<m;j++)
		{
			if (d[j]>0) /* place coefficients in the diagonal block */
			{
			d[j]++;
			w=int_matmatchrow(Ik.combs,Ik.card,Ik.n,d);
			t=d2diagblock[w];
			d[j]--;
			{
				mpz_t L_times_mi;
				mpz_init(L_times_mi);
				mpz_mul_ui(L_times_mi, L[j][r-1], mi[j]+d[j]);
				mpq_mspset_z(ls->B2r,lastB2r-r,w*r,L_times_mi,maxL[r-1]);
				mpz_clear(L_times_mi);
			}
			for (s=1;s<=r-1;s++)
			{
				mpz_t L_times_mi_r, neg_L_times_mi_s;
				mpz_init(L_times_mi_r);
				mpz_init(neg_L_times_mi_s);
				mpz_mul_ui(L_times_mi_r, L[j][r-1], mi[j]+d[j]);
				mpq_mspset_z(ls->B2r,lastB2r-r+s,w*r+s,L_times_mi_r,maxL[r-1]);
				mpz_mul_ui(neg_L_times_mi_s, L[j][s-1], mi[j]+d[j]);
				mpz_neg(neg_L_times_mi_s, neg_L_times_mi_s);
				mpq_t divisor;
				mpq_init(divisor);
				mpq_set_ui(divisor, maxL[s-1], 1);
				mpq_set_z(ls->C[t].diag[last[t]+s-1][d2relpos[w]*r+s], neg_L_times_mi_s);
				mpq_div(ls->C[t].diag[last[t]+s-1][d2relpos[w]*r+s], ls->C[t].diag[last[t]+s-1][d2relpos[w]*r+s], divisor);
				mpq_clear(divisor);
				mpz_clear(L_times_mi_r);
				mpz_clear(neg_L_times_mi_s);
			}
			}
		}
		for (j=0;j<m;j++)
		{
			if (d[j]==0) /* place coefficients in the non-diagonal block */
			{
			d[j]++;
			w=int_matmatchrow(Ik.combs,Ik.card,Ik.n,d);
			d[j]--;
			{
				mpz_t L_times_mi;
				mpz_init(L_times_mi);
				mpz_mul_ui(L_times_mi, L[j][r-1], mi[j]+d[j]);
				mpq_mspset_z(ls->B2r,lastB2r-r,w*r,L_times_mi,maxL[r-1]);
				mpz_clear(L_times_mi);
			}
			for (s=1;s<=r-1;s++)
			{
				mpz_t L_times_mi_r, neg_L_times_mi_s;
				mpz_init(L_times_mi_r);
				mpz_init(neg_L_times_mi_s);
				mpz_mul_ui(L_times_mi_r, L[j][r-1], mi[j]+d[j]);
				mpq_mspset_z(ls->B2r,lastB2r-r+s,w*r+s,L_times_mi_r,maxL[r-1]);
				mpz_mul_ui(neg_L_times_mi_s, L[j][s-1], mi[j]+d[j]);
				mpz_neg(neg_L_times_mi_s, neg_L_times_mi_s);
				mpq_mspset_z(ls->C[t].nondiag,last[t]+s-1,w*r+s,neg_L_times_mi_s,maxL[s-1]);
				mpz_clear(L_times_mi_r);
				mpz_clear(neg_L_times_mi_s);
			}
			
			}
		}
		for (s=1;s<=r-1;s++)
		{
			mpq_mspset_si(ls->A12,pcpos[lastpcpos],i*r,n[s-1],1);
			{
				mpz_t neg_z;
				mpz_init(neg_z);
				mpz_neg(neg_z, Z[s-1]);
				mpq_mspset_z(ls->A12,pcpos[lastpcpos++],i*r+s,neg_z,maxL[s-1]);
				mpz_clear(neg_z);
			}
		}
		w++;
		last[t]+=r-1;
	}
	t=0;
	for (h=1;h<=ls->H;h++)
	{
		for (k=0;k<nck(m,h);k++)
		{
//		if(DEBUG) 
		//printf("%d/%d\n ",k,nck(m,h)-1);
		fflush(stdout);
//		if(DEBUG) mpq_matprint(ls->C[t+k].diag,nck(r-1,r-h)*r,nck(r-1,r-h)*r);
//		if(DEBUG) mpq_mspprint(ls->C[t+k].nondiag);
		t0=CPUTIME;
		if (show_progress) {
			ls->C[t+k].lu_indices=(int*)mpq_ludcmp_progress(ls->C[t+k].diag,nck(r-1,r-h)*r, setup_start, n, r, show_progress);
		} else {
			ls->C[t+k].lu_indices=(int*)mpq_ludcmp(ls->C[t+k].diag,nck(r-1,r-h)*r);
		}
		if (ls->C[t+k].lu_indices == NULL) {
			// Singular matrix encountered during LU decomposition
			// Free allocated memory before returning
			int ii, jj;
			for (ii = 0; ii < t+k+1; ii++) {
				if (ls->C[ii].diag) {
					for (jj = 0; jj < nck(r-1,r-h)*r; jj++) {
						free(ls->C[ii].diag[jj]);
					}
					free(ls->C[ii].diag);
				}
				if (ls->C[ii].nondiag) {
					if (ls->C[ii].nondiag->coeff) {
						int jj;
						for (jj = 0; jj < ls->C[ii].nondiag->nnz; jj++) {
							mpq_clear(ls->C[ii].nondiag->coeff[jj]);
						}
						free(ls->C[ii].nondiag->coeff);
					}
					if (ls->C[ii].nondiag->pos_row)
						free(ls->C[ii].nondiag->pos_row);
					if (ls->C[ii].nondiag->pos_col)
						free(ls->C[ii].nondiag->pos_col);
					free(ls->C[ii].nondiag);
				}
			}
			free(ls->C);
			free(ls);
			free(Ik.combs);	
			free(I.combs);	
			return NULL;
		}
//		if(DEBUG) mpq_matprint(ls->C[t+k].diag,nck(r-1,r-h)*r,nck(r-1,r-h)*r);
		t1=CPUTIME;
		}
		t+=nck(m,h);
	}
	free(Ik.combs);	
	free(I.combs);	
return ls;
}

