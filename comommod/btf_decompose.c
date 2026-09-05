#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"

/* Transliteration of comom/btf_decompose.c into Z/p.  The block structure
 * is combinatorial, derived from the non-zero pattern of the Dn
 * combinations, so it is identical in every modular image; only the
 * "is this entry non-zero" test used to detect coupling changes meaning,
 * and an entry that vanishes mod p merely drops an off-diagonal block
 * that would have contributed nothing anyway. */

static int pattern_bitmask(int* comb, int r)
{
    int mask = 0;
    for (int s = 0; s < r - 1; s++)
        if (comb[s] > 0)
            mask |= (1 << s);
    return mask;
}

static int popcount(int mask)
{
    int count = 0;
    while (mask) { count += mask & 1; mask >>= 1; }
    return count;
}

static int find_block_index(int* block_masks, int num_blocks, int mask)
{
    for (int b = 0; b < num_blocks; b++)
        if (block_masks[b] == mask)
            return b;
    return -1;
}

btf_info* btf_decompose(zp_mat_t A11, combsrep* Dn, int r, int M, int cardGk)
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

        btf->blocks[0].diag = (zp_mat_t)zp_matzeros(cardGk, cardGk);
        for (int i = 0; i < cardGk; i++)
            memcpy(btf->blocks[0].diag[i], A11[i], (size_t)cardGk*sizeof(zp_t));

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

    /* Step 1: discover blocks from the Dn ordering */
    int total_combos = Dn->card;
    int* combo_mask = (int*)calloc(total_combos, sizeof(int));
    int* combo_filtered = (int*)calloc(total_combos, sizeof(int));

    for (int d = 0; d < total_combos; d++) {
        combo_mask[d] = pattern_bitmask(Dn->combs[d], r);
        combo_filtered[d] = (int_vecsubsum(Dn->combs[d], 0, r - 1) < M) ? 1 : 0;
    }

    int max_blocks = 1 << (r - 1);
    int* block_masks = (int*)calloc(max_blocks, sizeof(int));
    int num_blocks = 0;

    for (int d = 0; d < total_combos; d++) {
        int mask = combo_mask[d];
        if (find_block_index(block_masks, num_blocks, mask) < 0)
            block_masks[num_blocks++] = mask;
    }

    int* combo_block = (int*)calloc(total_combos, sizeof(int));
    for (int d = 0; d < total_combos; d++)
        combo_block[d] = find_block_index(block_masks, num_blocks, combo_mask[d]);

    /* Step 2: column permutation */
    btf->col_perm = (int*)calloc(cardGk, sizeof(int));
    int* block_ncols = (int*)calloc(num_blocks, sizeof(int));
    for (int d = 0; d < total_combos; d++)
        block_ncols[combo_block[d]] += M;

    int* block_col_offset = (int*)calloc(num_blocks + 1, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_col_offset[b + 1] = block_col_offset[b] + block_ncols[b];

    int* block_col_cursor = (int*)calloc(num_blocks, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_col_cursor[b] = block_col_offset[b];

    for (int d = 0; d < total_combos; d++) {
        int b = combo_block[d];
        for (int k = 0; k < M; k++)
            btf->col_perm[block_col_cursor[b]++] = d * M + k;
    }

    /* Step 3: row permutation, replicating setupls row generation order */
    int* block_nrows = (int*)calloc(num_blocks, sizeof(int));
    {
        for (int d = 0; d < total_combos; d++) {
            if (!combo_filtered[d]) continue;
            int b_src = combo_block[d];
            block_nrows[b_src] += M;
            for (int s = 0; s < r - 1; s++) {
                if (Dn->combs[d][s] > 0) {
                    block_nrows[b_src]++;
                } else {
                    int target_mask = combo_mask[d] | (1 << s);
                    int b_target = find_block_index(block_masks, num_blocks, target_mask);
                    if (b_target < 0) {
                        fprintf(stderr, "btf_decompose: cannot find block for mask 0x%x\n", target_mask);
                        block_nrows[b_src]++;
                    } else {
                        block_nrows[b_target]++;
                    }
                }
            }
        }
    }

    btf->row_perm = (int*)calloc(cardGk, sizeof(int));

    int* block_row_offset = (int*)calloc(num_blocks + 1, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_row_offset[b + 1] = block_row_offset[b] + block_nrows[b];

    int* block_row_cursor = (int*)calloc(num_blocks, sizeof(int));
    for (int b = 0; b < num_blocks; b++)
        block_row_cursor[b] = block_row_offset[b];

    {
        int row = 0;
        for (int d = 0; d < total_combos; d++) {
            if (!combo_filtered[d]) continue;
            int b_src = combo_block[d];

            for (int k = 0; k < M; k++) {
                btf->row_perm[block_row_cursor[b_src]++] = row;
                row++;
            }

            for (int s = 0; s < r - 1; s++) {
                int b_dest;
                if (Dn->combs[d][s] > 0) {
                    b_dest = b_src;
                } else {
                    int target_mask = combo_mask[d] | (1 << s);
                    b_dest = find_block_index(block_masks, num_blocks, target_mask);
                    if (b_dest < 0) b_dest = b_src;
                }
                btf->row_perm[block_row_cursor[b_dest]++] = row;
                row++;
            }
        }
    }

    /* Step 4: allocate and populate */
    btf->num_blocks = num_blocks;
    btf->blocks = (btf_block*)calloc(num_blocks, sizeof(btf_block));
    btf->block_starts = (int*)calloc(num_blocks + 1, sizeof(int));

    for (int b = 0; b < num_blocks; b++) {
        btf->blocks[b].size = block_ncols[b];
        btf->blocks[b].level = popcount(block_masks[b]);
        btf->blocks[b].lu_indices = NULL;
        btf->block_starts[b] = block_col_offset[b];
    }
    btf->block_starts[num_blocks] = cardGk;

    /* Step 5: diagonal sub-matrices */
    for (int b = 0; b < num_blocks; b++) {
        int sz = btf->blocks[b].size;
        if (sz == 0) { btf->blocks[b].diag = NULL; continue; }
        btf->blocks[b].diag = (zp_mat_t)zp_matzeros(sz, sz);
        int row_start = block_row_offset[b];
        int col_start = block_col_offset[b];
        for (int i = 0; i < sz && i < block_nrows[b]; i++) {
            int orig_row = btf->row_perm[row_start + i];
            for (int j = 0; j < sz; j++) {
                int orig_col = btf->col_perm[col_start + j];
                btf->blocks[b].diag[i][j] = A11[orig_row][orig_col];
            }
        }
    }

    /* Step 6: off-diagonal sub-matrices */
    int num_offdiag = 0;
    for (int b_row = 0; b_row < num_blocks; b_row++) {
        for (int b_col = b_row + 1; b_col < num_blocks; b_col++) {
            int mask_diff = block_masks[b_col] & ~block_masks[b_row];
            if (mask_diff == 0) continue;
            int has_nonzero = 0;
            int r_start = block_row_offset[b_row];
            int c_start = block_col_offset[b_col];
            for (int i = 0; i < block_nrows[b_row] && !has_nonzero; i++) {
                int orig_row = btf->row_perm[r_start + i];
                for (int j = 0; j < block_ncols[b_col] && !has_nonzero; j++) {
                    int orig_col = btf->col_perm[c_start + j];
                    if (A11[orig_row][orig_col] != 0) has_nonzero = 1;
                }
            }
            if (has_nonzero) num_offdiag++;
        }
    }

    btf->num_offdiag = num_offdiag;
    btf->offdiag = (num_offdiag > 0) ? (btf_offdiag*)calloc(num_offdiag, sizeof(btf_offdiag)) : NULL;

    int od_idx = 0;
    for (int b_row = 0; b_row < num_blocks; b_row++) {
        for (int b_col = b_row + 1; b_col < num_blocks; b_col++) {
            int mask_diff = block_masks[b_col] & ~block_masks[b_row];
            if (mask_diff == 0) continue;

            int r_start = block_row_offset[b_row];
            int c_start = block_col_offset[b_col];
            int nr = block_nrows[b_row];
            int nc = block_ncols[b_col];

            int has_nonzero = 0;
            for (int i = 0; i < nr && !has_nonzero; i++) {
                int orig_row = btf->row_perm[r_start + i];
                for (int j = 0; j < nc && !has_nonzero; j++) {
                    int orig_col = btf->col_perm[c_start + j];
                    if (A11[orig_row][orig_col] != 0) has_nonzero = 1;
                }
            }
            if (!has_nonzero) continue;

            btf->offdiag[od_idx].from_block = b_row;
            btf->offdiag[od_idx].to_block = b_col;
            btf->offdiag[od_idx].nrows = nr;
            btf->offdiag[od_idx].ncols = nc;
            btf->offdiag[od_idx].mat = (zp_mat_t)zp_matzeros(nr, nc);
            for (int i = 0; i < nr; i++) {
                int orig_row = btf->row_perm[r_start + i];
                for (int j = 0; j < nc; j++) {
                    int orig_col = btf->col_perm[c_start + j];
                    btf->offdiag[od_idx].mat[i][j] = A11[orig_row][orig_col];
                }
            }
            od_idx++;
        }
    }

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

    return btf;
}
