/*
 * linbox_solve.cpp — LinBox Dixon p-adic lifting wrapper for mpq_t systems.
 *
 * Converts the mpq_t linear system A*x = b into an integer system
 * by clearing denominators, solves via LinBox's Dixon solver, and
 * writes the rational solution back into b[].
 */

#include <cstdio>
#include <cstdlib>
#include <gmp.h>

/* Suppress unused-parameter warnings from LinBox/Givaro headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"

#include <givaro/zring.h>
#include <linbox/integer.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/solutions/solve.h>

#pragma GCC diagnostic pop

#include "linbox_solve.h"

extern "C" {

int linbox_solve_dense(mpq_mat_t A, mpq_vec_t b, int N)
{
    using Ring    = Givaro::ZRing<Givaro::Integer>;
    using Matrix  = LinBox::DenseMatrix<Ring>;
    using Vector  = LinBox::DenseVector<Ring>;

    Ring ZZ;

    /*
     * Convert the rational system A*x = b into an integer system.
     *
     * For each row i, compute the lcm of all denominators in that row
     * (including b[i]), then multiply through.  This gives an integer
     * matrix Aint and integer vector bint such that the rational solution
     * is the same.
     */
    Matrix Aint(ZZ, N, N);
    Vector bint(ZZ, N);

    mpz_t row_lcm, tmp_z;
    mpz_init(row_lcm);
    mpz_init(tmp_z);

    for (int i = 0; i < N; i++) {
        /* Compute lcm of all denominators in row i */
        mpz_set_ui(row_lcm, 1);
        for (int j = 0; j < N; j++) {
            mpq_canonicalize(A[i][j]);
            mpz_lcm(row_lcm, row_lcm, mpq_denref(A[i][j]));
        }
        mpq_canonicalize(b[i]);
        mpz_lcm(row_lcm, row_lcm, mpq_denref(b[i]));

        /* Multiply each entry by row_lcm to get integer entries */
        for (int j = 0; j < N; j++) {
            /* Aint[i][j] = A[i][j] * row_lcm = num(A[i][j]) * (row_lcm / den(A[i][j])) */
            mpz_divexact(tmp_z, row_lcm, mpq_denref(A[i][j]));
            mpz_mul(tmp_z, tmp_z, mpq_numref(A[i][j]));
            Givaro::Integer val;
            mpz_set(val.get_mpz(), tmp_z);
            Aint.setEntry(i, j, val);
        }

        /* bint[i] = b[i] * row_lcm */
        mpz_divexact(tmp_z, row_lcm, mpq_denref(b[i]));
        mpz_mul(tmp_z, tmp_z, mpq_numref(b[i]));
        Givaro::Integer val;
        mpz_set(val.get_mpz(), tmp_z);
        bint.setEntry(i, val);
    }

    mpz_clear(row_lcm);
    mpz_clear(tmp_z);

    /* Solve using Dixon p-adic lifting */
    Vector xNum(ZZ, N);
    Givaro::Integer xDen;

    try {
        LinBox::solve(xNum, xDen, Aint, bint,
                      LinBox::RingCategories::IntegerTag(),
                      LinBox::Method::Dixon());
    } catch (...) {
        return -1;
    }

    if (ZZ.isZero(xDen)) {
        return -1;
    }

    /* Write solution back into b[] as mpq_t values: b[i] = xNum[i] / xDen */
    for (int i = 0; i < N; i++) {
        mpz_set(mpq_numref(b[i]), xNum.getEntry(i).get_mpz_const());
        mpz_set(mpq_denref(b[i]), xDen.get_mpz_const());
        mpq_canonicalize(b[i]);
    }

    return 0;
}

} /* extern "C" */
