#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <fpla.h>
#include "momf.h"

/* Block back-substitution over the upper block triangular form of the
 * coefficient matrix.  Identical in structure to mom/blocksolve.c; the
 * only addition is the optional Wilkinson/Moler refinement of each
 * diagonal-block solve, which needs the pristine block kept in .orig. */
fp_vec_t blocksolve(LS* ls, fp_vec_t b)
{
	int h,k,t,i,w;

	fp_vec_t g=(fp_vec_t)fp_vec(nck(ls->m+ls->r-1,ls->r)*ls->r,0,1);
	t=1;
	w=0;
	for (h=ls->H;h>=1;h--)
	{
		int n=nck(ls->r-1,ls->r-h)*ls->r;
		fp_vec_t bg=(fp_vec_t)fp_vec(n,0,1);
		fp_vec_t xg=(fp_vec_t)fp_vec(n,0,1);
		for (k=nck(ls->m,h)-1;k>=0;k--)
		{
		if(DEBUG)printf("blocksolve: processing block %d/%lld with %d non-zeros \n",k+1,nck(ls->m,h),h);
		w+=n;
		fp_mspvecmul(bg,ls->C[ls->numdiagblocks-t].nondiag,g);
		for (i=n-1;i>=0;i--)
		{
			AORSCTR++;
			fp_sub(bg[i],b[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		}

		REFINE_SOLVES++;
		if (REFINE_ITERS > 0)
		{
			int sw = fp_lurefine(ls->C[ls->numdiagblocks-t].orig,
			                     ls->C[ls->numdiagblocks-t].diag,
			                     ls->C[ls->numdiagblocks-t].lu_indices,
			                     bg, xg, n, REFINE_ITERS, REFINE_TOL, NULL);
			if (sw < 0)
			{
				for (i=n;i>=0;i--) { fp_clear(bg[i]); fp_clear(xg[i]); }
				free(bg); free(xg);
				for (i=0;i<nck(ls->m+ls->r-1,ls->r)*ls->r;i++) fp_clear(g[i]);
				free(g);
				return NULL;
			}
			REFINE_SWEEPS += sw;
			fp_vecdup(bg, xg, n);
		}
		else if (fp_lubksb(ls->C[ls->numdiagblocks-t].diag,bg,n,ls->C[ls->numdiagblocks-t].lu_indices) < 0)
		{
			for (i=n;i>=0;i--) { fp_clear(bg[i]); fp_clear(xg[i]); }
			free(bg); free(xg);
			for (i=0;i<nck(ls->m+ls->r-1,ls->r)*ls->r;i++) fp_clear(g[i]);
			free(g);
			return NULL;
		}

		for (i=n-1;i>=0;i--)
			fp_set(g[nck(ls->m+ls->r-1,ls->r)*ls->r-w+i],bg[i]);
		t++;
		}

		for (i=n;i>=0;i--) { fp_clear(bg[i]); fp_clear(xg[i]); }
		free(bg); free(xg);
	}
	return g;
}
