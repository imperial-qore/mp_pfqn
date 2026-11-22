#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmpla.h"

/* set the denominator of an integer-valued q to den and keep q canonical
 * (GMP requires den > 0 and gcd(num,den)=1; mpz_set_ui turned a negative
 * den into a huge positive one). */
static void mpq_setden_checked(mpq_t q, int den)
{
	if (den == 0) { fprintf(stderr, "gmpla: zero denominator\n"); exit(1); }
	mpz_set_si(mpq_denref(q), den);
	mpq_canonicalize(q);
}

void mpq_mspset_z(mpq_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden)
{
	mpq_t* elem=(mpq_t*)mpq_mspget(msp,row,col);	
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;	
		msp->pos_col[msp->lastnnz]=col;	
		mpq_set_z(msp->coeff[msp->lastnnz], coeffnum);
		mpq_setden_checked(msp->coeff[msp->lastnnz], coeffden);
	}
	else
	{
		mpq_set_z(*elem, coeffnum);
		mpq_setden_checked(*elem, coeffden);
	}
}