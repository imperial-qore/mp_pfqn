#include "gmp.h"
#include "gmpla.h"

void mpq_mspset_z(mpq_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden)
{
	mpq_t* elem=(mpq_t*)mpq_mspget(msp,row,col);	
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;	
		msp->pos_col[msp->lastnnz]=col;	
		mpq_set_z(msp->coeff[msp->lastnnz], coeffnum);
		mpz_set_ui(mpq_denref(msp->coeff[msp->lastnnz]), coeffden);
	}
	else
	{
		mpq_set_z(*elem, coeffnum);
		mpz_set_ui(mpq_denref(*elem), coeffden);
	}
}