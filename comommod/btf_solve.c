#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"

/* btf_factorize - LU decompose each diagonal block in Z/p.
 *
 * Returns 0 on success, -1 if some block is singular modulo p.  Over a
 * finite field that is either a genuine singularity or an unlucky prime;
 * the caller falls back to a direct LU of A11 and, failing that, discards
 * the prime. */
int btf_factorize(btf_info* btf)
{
    if (!btf) return -1;

    for (int b = 0; b < btf->num_blocks; b++) {
        int sz = btf->blocks[b].size;
        if (sz == 0) continue;

        btf->blocks[b].lu_indices = zp_ludcmp(btf->blocks[b].diag, sz);
        if (btf->blocks[b].lu_indices == NULL)
            return -1;
    }
    return 0;
}

/* btf_solve - backward block substitution over the permuted UBT system.
 * The solution overwrites b in place.  Returns 0 on success, -1 on error. */
int btf_solve(btf_info* btf, zp_vec_t b, int cardGk)
{
    if (!btf) return -1;

    int n = btf->total_size;
    if (n != cardGk) {
        fprintf(stderr, "btf_solve: size mismatch (btf=%d, cardGk=%d)\n", n, cardGk);
        return -1;
    }

    zp_vec_t b_perm = (zp_vec_t)zp_vec(n, 0, 1);
    for (int i = 0; i < n; i++)
        b_perm[i] = b[btf->row_perm[i]];

    for (int blk = btf->num_blocks - 1; blk >= 0; blk--) {
        int sz = btf->blocks[blk].size;
        if (sz == 0) continue;
        int start = btf->block_starts[blk];

        /* subtract the contributions of the already-solved higher blocks */
        for (int od = 0; od < btf->num_offdiag; od++) {
            if (btf->offdiag[od].from_block != blk) continue;

            int j_blk = btf->offdiag[od].to_block;
            int j_start = btf->block_starts[j_blk];
            int j_sz = btf->blocks[j_blk].size;
            zp_mat_t mat = btf->offdiag[od].mat;
            int nr = btf->offdiag[od].nrows;

            for (int i = 0; i < nr; i++)
                for (int j = 0; j < j_sz; j++)
                    if (mat[i][j] != 0 && b_perm[j_start + j] != 0)
                        b_perm[start + i] = zp_subv(b_perm[start + i],
                                                    zp_mulv(mat[i][j], b_perm[j_start + j]));
        }

        /* solve the diagonal block in place */
        if (zp_lubksb(btf->blocks[blk].diag, &b_perm[start], sz,
                      btf->blocks[blk].lu_indices) < 0) {
            zp_vecfree(b_perm);
            return -1;
        }
    }

    for (int i = 0; i < n; i++)
        b[btf->col_perm[i]] = b_perm[i];

    zp_vecfree(b_perm);
    return 0;
}

void btf_free(btf_info* btf)
{
    if (!btf) return;

    if (btf->blocks) {
        for (int b = 0; b < btf->num_blocks; b++) {
            if (btf->blocks[b].diag)
                zp_matfree(btf->blocks[b].diag, btf->blocks[b].size);
            free(btf->blocks[b].lu_indices);
        }
        free(btf->blocks);
    }

    if (btf->offdiag) {
        for (int od = 0; od < btf->num_offdiag; od++)
            if (btf->offdiag[od].mat)
                zp_matfree(btf->offdiag[od].mat, btf->offdiag[od].nrows);
        free(btf->offdiag);
    }

    free(btf->row_perm);
    free(btf->col_perm);
    free(btf->block_starts);
    free(btf);
}
