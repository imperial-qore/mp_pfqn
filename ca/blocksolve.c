#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "mom.h"

mpq_vec_t blocksolve(LS* ls, mpq_vec_t b)
{
	int h,k,t,i,w;
#ifdef VERBOSE_LEVEL_2
	printf("entering blocksolve. \n");
#endif
	
	mpq_vec_t g=(mpq_vec_t)mpq_vec(nck(ls->m+ls->r-1,ls->r)*ls->r,0,1);
	t=1;
	w=0;
	for (h=ls->H;h>=1;h--)
	{
		int n=nck(ls->r-1,ls->r-h)*ls->r;
		mpq_vec_t bg=(mpq_vec_t)mpq_vec(n,0,1);
		for (k=nck(ls->m,h)-1;k>=0;k--)
		{
		if(DEBUG)printf("blocksolve: processing block %d/%lld with %d non-zeros \n",k+1,nck(ls->m,h),h);
		w+=n;
		mpq_mspvecmul(bg,ls->C[ls->numdiagblocks-t].nondiag,g);	
		for (i=n-1;i>=0;i--)
		{
			AORSCTR++;
			mpq_sub(bg[i],b[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		}
		mpq_lubksb(ls->C[ls->numdiagblocks-t].diag,bg,n,ls->C[ls->numdiagblocks-t].lu_indices);
		
		for (i=n-1;i>=0;i--)
			mpq_set(g[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		t++;
		}
		
		for (i=n-1;i>=0;i--)
			mpq_clear(bg[i]);
	}
	return g;
}
