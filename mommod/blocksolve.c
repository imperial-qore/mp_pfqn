#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <zpla.h>
#include "mommod.h"

/* Block back-substitution over the upper block triangular form of the
 * coefficient matrix, carried out entirely in Z/p.  Identical in
 * structure to mom/blocksolve.c. */
zp_vec_t blocksolve(LS* ls, zp_vec_t b)
{
	int h,k,t,i,w;

	zp_vec_t g=(zp_vec_t)zp_vec(nck(ls->m+ls->r-1,ls->r)*ls->r,0,1);
	t=1;
	w=0;
	for (h=ls->H;h>=1;h--)
	{
		int n=nck(ls->r-1,ls->r-h)*ls->r;
		zp_vec_t bg=(zp_vec_t)zp_vec(n,0,1);
		for (k=nck(ls->m,h)-1;k>=0;k--)
		{
		if(DEBUG)printf("blocksolve: processing block %d/%lld with %d non-zeros \n",k+1,nck(ls->m,h),h);
		w+=n;
		zp_mspvecmul(bg,ls->C[ls->numdiagblocks-t].nondiag,g);
		for (i=n-1;i>=0;i--)
		{
			AORSCTR++;
			zp_sub(bg[i],b[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		}

		if (zp_lubksb(ls->C[ls->numdiagblocks-t].diag,bg,n,ls->C[ls->numdiagblocks-t].lu_indices) < 0)
		{
			zp_vecfree(bg); zp_vecfree(g);
			return NULL;
		}

		for (i=n-1;i>=0;i--)
			zp_set(g[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		t++;
		}

		zp_vecfree(bg);
	}
	return g;
}
