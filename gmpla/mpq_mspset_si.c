#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmpla.h"

/* mpq_set_si takes an unsigned denominator: a negative coeffden became a
 * huge value and the result was not canonical. */
static void mpq_set_si_checked(mpq_t q, long num, long den)
{
	if (den == 0) { fprintf(stderr, "gmpla: zero denominator\n"); exit(1); }
	if (den < 0) { num = -num; den = -den; }
	mpq_set_si(q, num, (unsigned long) den);
	mpq_canonicalize(q);
}

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
		mpq_set_si_checked(msp->coeff[msp->lastnnz],coeffnum,coeffden);
	}
	else
	{
		mpq_set_si_checked(*elem,coeffnum,coeffden);
	}
}
