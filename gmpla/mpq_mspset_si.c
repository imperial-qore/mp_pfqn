#include "gmp.h"
#include "gmpla.h"

void mpq_mspset_si(mpq_msp_t msp, int row, int col, int coeffnum, int coeffden)
{
	mpq_t* elem=(mpq_t*)mpq_mspget(msp,row,col);	
	if (elem == NULL)
	{
		msp->lastnnz++;
/*		if (row+1 > msp->rows)
		{
			printf("mpq_mspsetsi: error - row index outside allowed range\n");
			return;
		}
		if (col+1 > msp->cols)
		{
			printf("mpq_mspsetsi: error - col index outside allowed range\n");
			mpq_mspprint(msp);
			exit(1);
		}
		if (msp->lastnnz+1>msp->nnz)
		{
			printf("mpq_mspsetsi: error - maximum number of non-zeros exceeded\n");
			return;
		}
*/		msp->pos_row[msp->lastnnz]=row;	
		msp->pos_col[msp->lastnnz]=col;	
		mpq_set_si(msp->coeff[msp->lastnnz],coeffnum,coeffden);
	}
	else
	{
		mpq_set_si(*elem,coeffnum,coeffden);
	}
}
