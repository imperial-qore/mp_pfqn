#include "gmpla.h"

/* Matrix transpose multiply: rop = A' * B
   A is nrows x ncols, B is nrows x ncols, rop is ncols x ncols */
void mpq_mattransmul(mpq_mat_t rop, mpq_mat_t A, mpq_mat_t B, int nrows, int ncols)
{
	int i, j, k;
	mpq_t t;
	mpq_init(t);

	for (i = 0; i < ncols; i++) {
		for (j = 0; j < ncols; j++) {
			mpq_set_si(rop[i][j], 0, 1);
			for (k = 0; k < nrows; k++) {
				if (mpq_sgn(A[k][i]) != 0 && mpq_sgn(B[k][j]) != 0) {
					mpq_mul(t, A[k][i], B[k][j]);
					mpq_add(rop[i][j], rop[i][j], t);
				}
			}
		}
	}

	mpq_clear(t);
}
