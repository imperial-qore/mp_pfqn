#include "gmpla.h"

/* Dense matrix-vector multiply: rop = A * v, where A is N x N */
void mpq_matvecmul(mpq_vec_t rop, mpq_mat_t A, mpq_vec_t v, int N)
{
	int i, j;
	mpq_t t;
	mpq_init(t);

	for (i = 0; i < N; i++) {
		mpq_set_si(rop[i], 0, 1);
		for (j = 0; j < N; j++) {
			if (mpq_sgn(A[i][j]) != 0 && mpq_sgn(v[j]) != 0) {
				mpq_mul(t, A[i][j], v[j]);
				mpq_add(rop[i], rop[i], t);
			}
		}
	}

	mpq_clear(t);
}
