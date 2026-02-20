#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"

/**
 * Compute bitmask of non-zero positions in comb[0..r-2].
 * This determines the block pattern for a combination.
 */
static int pattern_bitmask(int* comb, int r)
{
    int mask = 0;
    for (int s = 0; s < r - 1; s++) {
        if (comb[s] > 0)
            mask |= (1 << s);
    }
    return mask;
}

/**
 * Population count (number of set bits) = level of a block.
 */
static int popcount(int mask)
{
    int count = 0;
    while (mask) {
        count += mask & 1;
        mask >>= 1;
    }
    return count;
}

/**
 * Find the index of a bitmask in the block_masks array.
 * Returns -1 if not found.
 */
static int find_block_index(int* block_masks, int num_blocks, int mask)
{
    for (int b = 0; b < num_blocks; b++) {
        if (block_masks[b] == mask)
            return b;
    }
    return -1;
}

/**
 * btf_decompose - Discover block structure and extract sub-matrices from A11.
 *
 * A11 has upper block-triangular structure when rows/columns are permuted
 * according to the non-zero pattern of the Dn combinations. Combinations
 * sorted by sortbynnzpos naturally group into blocks by pattern bitmask,
 * ordered by increasing popcount (level). Coupling goes from lower to
 * higher levels = UBT.
 *
 * PC rows for class s NOT in the source combo's pattern are reassigned
 * to the block with pattern P|{s}, making all diagonal blocks square.
 */
btf_info* btf_decompose(mpq_mat_t A11, combsrep* Dn, int r, int M, int cardGk)
{
    btf_info* btf = (btf_info*)calloc(1, sizeof(btf_info));
    if (!btf) return NULL;
    btf->total_size = cardGk;

    /* --- Special case: r == 1 (single class, no PC rows) --- */
    if (r == 1) {
        btf->num_blocks = 1;
        btf->blocks = (btf_block*)calloc(1, sizeof(btf_block));
        btf->block_starts = (int*)calloc(2, sizeof(int));
        btf->block_starts[0] = 0;
        btf->block_starts[1] = cardGk;
        btf->blocks[0].size = cardGk;
        btf->blocks[0].level = 0;
        btf->blocks[0].lu_indices = NULL;

        /* Diagonal block = full A11 (copy) */
        btf->blocks[0].diag = (mpq_mat_t)mpq_matzeros(cardGk, cardGk);
        for (int i = 0; i < cardGk; i++)
            for (int j = 0; j < cardGk; j++)
                mpq_set(btf->blocks[0].diag[i][j], A11[i][j]);

        /* Identity permutations */
        btf->row_perm = (int*)calloc(cardGk, sizeof(int));
        btf->col_perm = (int*)calloc(cardGk, sizeof(int));
        for (int i = 0; i < cardGk; i++) {
            btf->row_perm[i] = i;
            btf->col_perm[i] = i;
        }

        btf->num_offdiag = 0;
        btf->offdiag = NULL;
        return btf;
    }

    /* ========================================================
     * Step 1: Discover blocks from Dn ordering
     * ======================================================== */

    /* Determine which combos are filtered and compute their bitmask.
     * Also assign ALL combos (filtered or not) to blocks by bitmask. */
    int total_combos = Dn->card;
    int* combo_mask = (int*)calloc(total_combos, sizeof(int));
    int* combo_filtered = (int*)calloc(total_combos, sizeof(int)); /* 1 if filtered */

    for (int d = 0; d < total_combos; d++) {
        combo_mask[d] = pattern_bitmask(Dn->combs[d], r);
        /* Filter: same as setupls.c line 40 */
        combo_filtered[d] = (int_vecsubsum(Dn->combs[d], 0, r - 1) < M) ? 1 : 0;
    }

    /* Collect distinct bitmasks in Dn order (they should already be
     * grouped by pattern since Dn is sorted by sortbynnzpos).
     * The maximum number of distinct patterns is 2^(r-1). */
    int max_blocks = 1 << (r - 1);
    int* block_masks = (int*)calloc(max_blocks, sizeof(int));
    int num_blocks = 0;

    for (int d = 0; d < total_combos; d++) {
        int mask = combo_mask[d];
        if (find_block_index(block_masks, num_blocks, mask) < 0) {
            block_masks[num_blocks++] = mask;
        }
    }

    /* Map each combo to its block index */
    int* combo_block = (int*)calloc(total_combos, sizeof(int));
    for (int d = 0; d < total_combos; d++) {
        combo_block[d] = find_block_index(block_masks, num_blocks, combo_mask[d]);
    }

    /* ========================================================
     * Step 2: Build column permutation
     *
     * Columns are ordered block by block. Within each block,
     * combos appear in their Dn order, each contributing M columns.
     * ======================================================== */

    btf->col_perm = (int*)calloc(cardGk, sizeof(int));
    int* block_ncols = (int*)calloc(num_blocks, sizeof(int));

    /* Count columns per block */
    for (int d = 0; d < total_combos; d++) {
        block_ncols[combo_block[d]] += M;
    }

    /* Compute column offsets per block */
    int* block_col_offset = (int*)calloc(num_blocks + 1, sizeof(int));
    for (int b = 0; b < num_blocks; b++) {
        block_col_offset[b + 1] = block_col_offset[b] + block_ncols[b];
    }

    /* Fill col_perm: for each block, for each combo in that block (in Dn order),
     * map M columns. We need a running counter per block. */
    int* block_col_cursor = (int*)calloc(num_blocks, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_col_cursor[b] = block_col_offset[b];

    for (int d = 0; d < total_combos; d++) {
        int b = combo_block[d];
        /* The M original columns for combo d are at positions d*M + 0 .. d*M + M-1
         * (because hash(Dn,comb,k) = d*M + (k-1) + 1, so 0-indexed = d*M + k-1) */
        for (int k = 0; k < M; k++) {
            btf->col_perm[block_col_cursor[b]++] = d * M + k;
        }
    }

    /* ========================================================
     * Step 3: Build row permutation
     *
     * Replicate setupls.c row generation order exactly:
     * For each filtered combo d:
     *   M CE rows → stay in source block
     *   r-1 PC rows:
     *     s in P (comb[s]>0) → stay in source block
     *     s not in P → move to block P|{s}
     * ======================================================== */

    /* First pass: count rows per block */
    int* block_nrows = (int*)calloc(num_blocks, sizeof(int));
    {
        int row = 0;
        for (int d = 0; d < total_combos; d++) {
            if (!combo_filtered[d]) continue;
            int b_src = combo_block[d];

            /* M CE rows → source block */
            block_nrows[b_src] += M;
            row += M;

            /* r-1 PC rows */
            for (int s = 0; s < r - 1; s++) {
                if (Dn->combs[d][s] > 0) {
                    /* s in P → stays */
                    block_nrows[b_src]++;
                } else {
                    /* s not in P → goes to block P|{s} */
                    int target_mask = combo_mask[d] | (1 << s);
                    int b_target = find_block_index(block_masks, num_blocks, target_mask);
                    if (b_target < 0) {
                        fprintf(stderr, "btf_decompose: cannot find block for mask 0x%x\n", target_mask);
                        /* fallback */
                        block_nrows[b_src]++;
                    } else {
                        block_nrows[b_target]++;
                    }
                }
                row++;
            }
        }
    }

    /* Build row_perm by collecting rows block-by-block */
    btf->row_perm = (int*)calloc(cardGk, sizeof(int));

    /* Compute row offsets per block (same as block_starts) */
    int* block_row_offset = (int*)calloc(num_blocks + 1, sizeof(int));
    for (int b = 0; b < num_blocks; b++) {
        block_row_offset[b + 1] = block_row_offset[b] + block_nrows[b];
    }

    /* Fill row_perm using running cursors per block */
    int* block_row_cursor = (int*)calloc(num_blocks, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_row_cursor[b] = block_row_offset[b];

    /* Temporary array: for each original row, which block it goes to */
    int* row_to_block = (int*)calloc(cardGk, sizeof(int));

    {
        int row = 0;
        for (int d = 0; d < total_combos; d++) {
            if (!combo_filtered[d]) continue;
            int b_src = combo_block[d];

            /* M CE rows → source block */
            for (int k = 0; k < M; k++) {
                row_to_block[row] = b_src;
                btf->row_perm[block_row_cursor[b_src]++] = row;
                row++;
            }

            /* r-1 PC rows */
            for (int s = 0; s < r - 1; s++) {
                int b_dest;
                if (Dn->combs[d][s] > 0) {
                    b_dest = b_src;
                } else {
                    int target_mask = combo_mask[d] | (1 << s);
                    b_dest = find_block_index(block_masks, num_blocks, target_mask);
                    if (b_dest < 0) b_dest = b_src;
                }
                row_to_block[row] = b_dest;
                btf->row_perm[block_row_cursor[b_dest]++] = row;
                row++;
            }
        }
    }

    /* ========================================================
     * Step 4: Allocate and populate btf_info
     * ======================================================== */

    btf->num_blocks = num_blocks;
    btf->blocks = (btf_block*)calloc(num_blocks, sizeof(btf_block));
    btf->block_starts = (int*)calloc(num_blocks + 1, sizeof(int));

    for (int b = 0; b < num_blocks; b++) {
        /* Diagonal blocks must be square: nrows == ncols */
        if (block_nrows[b] != block_ncols[b]) {
            fprintf(stderr, "btf_decompose: block %d not square (rows=%d, cols=%d, mask=0x%x)\n",
                    b, block_nrows[b], block_ncols[b], block_masks[b]);
            /* Continue anyway; the solver will likely detect singularity */
        }
        btf->blocks[b].size = block_ncols[b]; /* use column count as block size */
        btf->blocks[b].level = popcount(block_masks[b]);
        btf->blocks[b].lu_indices = NULL;
        btf->block_starts[b] = block_col_offset[b];
    }
    btf->block_starts[num_blocks] = cardGk;

    /* ========================================================
     * Step 5: Extract diagonal sub-matrices
     * ======================================================== */

    for (int b = 0; b < num_blocks; b++) {
        int sz = btf->blocks[b].size;
        if (sz == 0) {
            btf->blocks[b].diag = NULL;
            continue;
        }
        btf->blocks[b].diag = (mpq_mat_t)mpq_matzeros(sz, sz);
        int row_start = block_row_offset[b];
        int col_start = block_col_offset[b];
        for (int i = 0; i < sz && i < block_nrows[b]; i++) {
            for (int j = 0; j < sz; j++) {
                int orig_row = btf->row_perm[row_start + i];
                int orig_col = btf->col_perm[col_start + j];
                mpq_set(btf->blocks[b].diag[i][j], A11[orig_row][orig_col]);
            }
        }
    }

    /* ========================================================
     * Step 6: Extract off-diagonal sub-matrices
     *
     * For each pair (b_row, b_col) where b_row < b_col and there
     * exists coupling from b_row's rows to b_col's columns:
     * pattern(b_col) = pattern(b_row) | {s} for some s.
     * ======================================================== */

    /* First count the number of off-diagonal blocks */
    int num_offdiag = 0;
    for (int b_row = 0; b_row < num_blocks; b_row++) {
        for (int b_col = b_row + 1; b_col < num_blocks; b_col++) {
            /* Check if b_col's mask is a superset of b_row's mask
             * differing by exactly one bit, OR any coupling exists */
            int mask_diff = block_masks[b_col] & ~block_masks[b_row];
            if (mask_diff == 0) continue; /* b_col not a superset of b_row */
            /* Any non-zero coupling? Check */
            int has_nonzero = 0;
            int r_start = block_row_offset[b_row];
            int c_start = block_col_offset[b_col];
            for (int i = 0; i < block_nrows[b_row] && !has_nonzero; i++) {
                for (int j = 0; j < block_ncols[b_col] && !has_nonzero; j++) {
                    int orig_row = btf->row_perm[r_start + i];
                    int orig_col = btf->col_perm[c_start + j];
                    if (mpq_sgn(A11[orig_row][orig_col]) != 0)
                        has_nonzero = 1;
                }
            }
            if (has_nonzero) num_offdiag++;
        }
    }

    btf->num_offdiag = num_offdiag;
    btf->offdiag = (num_offdiag > 0) ? (btf_offdiag*)calloc(num_offdiag, sizeof(btf_offdiag)) : NULL;

    /* Second pass: extract off-diagonal sub-matrices */
    int od_idx = 0;
    for (int b_row = 0; b_row < num_blocks; b_row++) {
        for (int b_col = b_row + 1; b_col < num_blocks; b_col++) {
            int mask_diff = block_masks[b_col] & ~block_masks[b_row];
            if (mask_diff == 0) continue;

            int r_start = block_row_offset[b_row];
            int c_start = block_col_offset[b_col];
            int nr = block_nrows[b_row];
            int nc = block_ncols[b_col];

            /* Check for any non-zero entry */
            int has_nonzero = 0;
            for (int i = 0; i < nr && !has_nonzero; i++) {
                for (int j = 0; j < nc && !has_nonzero; j++) {
                    int orig_row = btf->row_perm[r_start + i];
                    int orig_col = btf->col_perm[c_start + j];
                    if (mpq_sgn(A11[orig_row][orig_col]) != 0)
                        has_nonzero = 1;
                }
            }
            if (!has_nonzero) continue;

            btf->offdiag[od_idx].from_block = b_row;
            btf->offdiag[od_idx].to_block = b_col;
            btf->offdiag[od_idx].nrows = nr;
            btf->offdiag[od_idx].ncols = nc;
            btf->offdiag[od_idx].mat = (mpq_mat_t)mpq_matzeros(nr, nc);
            for (int i = 0; i < nr; i++) {
                for (int j = 0; j < nc; j++) {
                    int orig_row = btf->row_perm[r_start + i];
                    int orig_col = btf->col_perm[c_start + j];
                    mpq_set(btf->offdiag[od_idx].mat[i][j], A11[orig_row][orig_col]);
                }
            }
            od_idx++;
        }
    }

    /* Clean up temporaries */
    free(combo_mask);
    free(combo_filtered);
    free(combo_block);
    free(block_masks);
    free(block_ncols);
    free(block_nrows);
    free(block_col_offset);
    free(block_row_offset);
    free(block_col_cursor);
    free(block_row_cursor);
    free(row_to_block);

    return btf;
}
