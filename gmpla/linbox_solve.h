#ifndef LINBOX_SOLVE_H
#define LINBOX_SOLVE_H

#include <gmp.h>

/* mpq_t types from gmpla */
#ifndef MPQ_MAT
typedef mpq_t* mpq_vec_t;
typedef mpq_t** mpq_mat_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * linbox_solve_dense — Solve A*x = b using LinBox Dixon p-adic lifting.
 *
 * A is an N×N matrix of mpq_t (NOT modified).
 * b is an N-vector of mpq_t (overwritten with solution x on success).
 * Returns 0 on success, -1 on singular/failure.
 */
int linbox_solve_dense(mpq_mat_t A, mpq_vec_t b, int N);

#ifdef __cplusplus
}
#endif

#endif /* LINBOX_SOLVE_H */
