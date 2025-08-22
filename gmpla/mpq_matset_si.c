#include "gmp.h"
#include "gmpla.h"

void mpq_matset_si(mpq_mat_t mat, int row, int col, int coeffnum, int coeffden)
{
	mpq_set_si(mat[row][col],coeffnum,coeffden);
}
