#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"

/**
 * btf_factorize - LU decompose each diagonal block.
 *
 * Returns 0 on success, -1 if any block is singular.
 */
int btf_factorize(btf_info* btf)
{
    if (!btf) return -1;

    for (int b = 0; b < btf->num_blocks; b++) {
        int sz = btf->blocks[b].size;
        if (sz == 0) continue;

        btf->blocks[b].lu_indices = mpq_ludcmp(btf->blocks[b].diag, sz);
        if (btf->blocks[b].lu_indices == NULL)
            return -1; /* singular block */
    }
    return 0;
}

/**
 * btf_solve - Solve the permuted UBT system using backward block substitution.
 *
 * The system P_r * A11 * P_c^T has upper block-triangular structure.
 * We solve from the last (highest-level) block to the first:
 *
 *   1. Permute RHS:  b_perm[i] = b[row_perm[i]]
 *   2. For blk = num_blocks-1 downto 0:
 *        rhs = b_perm[block_starts[blk] .. +size-1]
 *        For each offdiag with from_block == blk:
 *            rhs -= offdiag.mat * x_perm[block_starts[to_block]..+size_to]
 *        Solve: mpq_lubksb(diag, rhs, size, lu_indices)
 *        Store rhs into x_perm
 *   3. Inverse column permutation:  b[col_perm[i]] = x_perm[i]
 *
 * The solution overwrites b in-place. Returns 0 on success, -1 on error.
 */
int btf_solve(btf_info* btf, mpq_vec_t b, int cardGk)
{
    if (!btf) return -1;

    int n = btf->total_size;
    if (n != cardGk) {
        fprintf(stderr, "btf_solve: size mismatch (btf=%d, cardGk=%d)\n", n, cardGk);
        return -1;
    }

    /* Step 1: Apply row permutation to RHS */
    mpq_vec_t b_perm = (mpq_vec_t)mpq_vec(n, 0, 1);
    for (int i = 0; i < n; i++)
        mpq_set(b_perm[i], b[btf->row_perm[i]]);

    /* Step 2: Backward block substitution */
    mpq_t tmp;
    mpq_init(tmp);

    for (int blk = btf->num_blocks - 1; blk >= 0; blk--) {
        int sz = btf->blocks[blk].size;
        if (sz == 0) continue;
        int start = btf->block_starts[blk];

        /* Subtract off-diagonal contributions from already-solved higher blocks */
        for (int od = 0; od < btf->num_offdiag; od++) {
            if (btf->offdiag[od].from_block != blk) continue;

            int j_blk = btf->offdiag[od].to_block;
            int j_start = btf->block_starts[j_blk];
            int j_sz = btf->blocks[j_blk].size;
            mpq_mat_t mat = btf->offdiag[od].mat;
            int nr = btf->offdiag[od].nrows;

            /* rhs[i] -= sum_j mat[i][j] * x_perm[j_start + j] */
            for (int i = 0; i < nr; i++) {
                for (int j = 0; j < j_sz; j++) {
                    if (mpq_sgn(mat[i][j]) != 0 && mpq_sgn(b_perm[j_start + j]) != 0) {
                        mpq_mul(tmp, mat[i][j], b_perm[j_start + j]);
                        mpq_sub(b_perm[start + i], b_perm[start + i], tmp);
                    }
                }
            }
        }

        /* Solve diagonal block */
        /* mpq_lubksb operates on a contiguous sub-vector; we need to
         * provide a pointer to the start of this block's portion. */
        mpq_vec_t block_rhs = &b_perm[start];

        /* mpq_lubksb expects the full LU matrix and the vector starting
         * at index 0. We need a temporary copy since lubksb uses 0-based
         * indexing into the vector. */
        mpq_vec_t rhs_tmp = (mpq_vec_t)mpq_vec(sz, 0, 1);
        for (int i = 0; i < sz; i++)
            mpq_set(rhs_tmp[i], block_rhs[i]);

        if (mpq_lubksb(btf->blocks[blk].diag, rhs_tmp, sz, btf->blocks[blk].lu_indices) < 0) {
            /* Singular block */
            for (int i = 0; i < sz; i++) mpq_clear(rhs_tmp[i]);
            free(rhs_tmp);
            for (int i = 0; i < n; i++) mpq_clear(b_perm[i]);
            free(b_perm);
            mpq_clear(tmp);
            return -1;
        }

        /* Copy solution back */
        for (int i = 0; i < sz; i++)
            mpq_set(block_rhs[i], rhs_tmp[i]);

        for (int i = 0; i < sz; i++) mpq_clear(rhs_tmp[i]);
        free(rhs_tmp);
    }

    /* Step 3: Apply inverse column permutation: b[col_perm[i]] = b_perm[i] */
    for (int i = 0; i < n; i++)
        mpq_set(b[btf->col_perm[i]], b_perm[i]);

    /* Clean up */
    for (int i = 0; i < n; i++) mpq_clear(b_perm[i]);
    free(b_perm);
    mpq_clear(tmp);

    return 0;
}

/**
 * btf_free - Free all memory associated with a btf_info structure.
 */
void btf_free(btf_info* btf)
{
    if (!btf) return;

    /* Free diagonal blocks */
    if (btf->blocks) {
        for (int b = 0; b < btf->num_blocks; b++) {
            int sz = btf->blocks[b].size;
            if (btf->blocks[b].diag) {
                for (int i = 0; i < sz; i++) {
                    if (btf->blocks[b].diag[i]) {
                        for (int j = 0; j < sz; j++)
                            mpq_clear(btf->blocks[b].diag[i][j]);
                        free(btf->blocks[b].diag[i]);
                    }
                }
                free(btf->blocks[b].diag);
            }
            if (btf->blocks[b].lu_indices)
                free(btf->blocks[b].lu_indices);
        }
        free(btf->blocks);
    }

    /* Free off-diagonal blocks */
    if (btf->offdiag) {
        for (int od = 0; od < btf->num_offdiag; od++) {
            int nr = btf->offdiag[od].nrows;
            int nc = btf->offdiag[od].ncols;
            if (btf->offdiag[od].mat) {
                for (int i = 0; i < nr; i++) {
                    if (btf->offdiag[od].mat[i]) {
                        for (int j = 0; j < nc; j++)
                            mpq_clear(btf->offdiag[od].mat[i][j]);
                        free(btf->offdiag[od].mat[i]);
                    }
                }
                free(btf->offdiag[od].mat);
            }
        }
        free(btf->offdiag);
    }

    /* Free permutations */
    if (btf->row_perm) free(btf->row_perm);
    if (btf->col_perm) free(btf->col_perm);
    if (btf->block_starts) free(btf->block_starts);

    free(btf);
}
