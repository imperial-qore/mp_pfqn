#include "gmp.h"
#include "gmpla.h"

mpq_t* mpq_mspget(mpq_msp_t msp, int row, int col)
{
	int k;
	for (k=0;k<=msp->lastnnz;k++)
		if (msp->pos_row[k]==row && msp->pos_col[k]==col) return &msp->coeff[k];
	return NULL;
}
