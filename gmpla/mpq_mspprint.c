#include <stdio.h>
#include <gmp.h>
#include "gmpla.h"

void mpq_mspprint(mpq_msp_t msp)
{
	int i,j;
	mpq_t* t;
	// find max row
	printf("max nnz=%d ",msp->nnz);
	printf("current nnz=%d\n",msp->lastnnz+1);
	for(i=0;i<msp->rows;i++)
	{
		for(j=0;j<msp->cols;j++)
		{
			t=(mpq_t*)mpq_mspget(msp,i,j);
			
			if (t != NULL) printf("%3g ",(double)mpq_get_d(*t));
			else printf("  0 ");
				
		}
		printf("\n");
	}
	printf("\n");
}
