#include "gmpla.h"

mpq_msp_t mpq_msp(int rows, int cols, int numnnz)
{
	mpq_msp_t msp=(mpq_msp_t)calloc(1,(sizeof(struct struct_mpq_msp)));
	msp->coeff=(mpq_vec_t) mpq_vec(numnnz,0,1);
	msp->pos_row=(int*) calloc((numnnz),sizeof(int));
	msp->pos_col=(int*) calloc((numnnz),sizeof(int));
	msp->nnz=numnnz;
	msp->lastnnz=-1;
	msp->rows=rows;
	msp->cols=cols;
	return msp;
}
