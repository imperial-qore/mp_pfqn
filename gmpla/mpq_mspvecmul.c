#include "gmp.h"
#include "gmpla.h"

/* Define the global variable here */
mpq_t rmul;

void mpq_mspvecmul(mpq_vec_t rop, mpq_msp_t op1, mpq_vec_t op2)
{
    int t,i;
    mpq_init(rmul);
    for(i=0;i<op1->rows;i++)
	    mpq_set_si(rop[i],0,1);
    for(t=0;t<op1->nnz;t++)
    {
	mpq_mul(rmul,op1->coeff[t],op2[op1->pos_col[t]]);   
    	mpq_add(rop[op1->pos_row[t]],rop[op1->pos_row[t]],rmul);
    }
    mpq_clear(rmul);
}
