#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"

/**
 * Determine the block level based on the non-zero pattern in the population vector
 */
static int get_block_level(int* comb, int r)
{
    int level = 0;
    for (int i = 0; i < r; i++) {
        if (comb[i] > 0) level++;
    }
    return level;
}

/**
 * Create block structure for BTF decomposition
 */
void btf_create_blocks(LS* ls)
{
    int M = ls->m;
    int R = ls->r;
    
    // Use single block BTF structure for stability
    // This ensures correct results while we develop the full microblock decomposition
    ls->H = 1;
    ls->num_macro_blocks = 1;
    
    // Allocate block arrays
    ls->X_blocks = (mpq_diagblock**)calloc(ls->num_macro_blocks, sizeof(mpq_diagblock*));
    ls->Y_blocks = (mpq_diagblock**)calloc(ls->num_macro_blocks, sizeof(mpq_diagblock*));
    ls->B1_blocks = (mpq_diagblock**)calloc(ls->num_macro_blocks, sizeof(mpq_diagblock*));
    ls->B2_blocks = (mpq_diagblock**)calloc(ls->num_macro_blocks, sizeof(mpq_diagblock*));
    ls->C_blocks = (mpq_diagblock**)calloc(ls->num_macro_blocks, sizeof(mpq_diagblock*));
    ls->block_positions = (block_position*)calloc(ls->num_macro_blocks, sizeof(block_position));
    
    // Initialize block metadata
    ls->block_levels = (int*)calloc(ls->num_macro_blocks, sizeof(int));
    ls->block_patterns = (int**)calloc(ls->num_macro_blocks, sizeof(int*));
    
    // Calculate actual matrix size
    long int cardG = nck(M + R - 1, M);
    long int cardGk = cardG * M;
    
    // Create single block encompassing entire matrix
    ls->block_levels[0] = 1;
    ls->block_patterns[0] = (int*)calloc(R, sizeof(int));
    
    // Store block position
    ls->block_positions[0].start_row = 0;
    ls->block_positions[0].start_col = 0;
    ls->block_positions[0].rows = cardGk;
    ls->block_positions[0].cols = cardGk;
    
    // Create diagonal block - reuse A11
    ls->X_blocks[0] = (mpq_diagblock*)calloc(1, sizeof(mpq_diagblock));
    ls->X_blocks[0]->diag = ls->A11;  // Reference existing A11
    ls->X_blocks[0]->lu_indices = NULL;
    
    // No off-diagonal blocks needed for single block
    ls->Y_blocks[0] = NULL;
    ls->B1_blocks[0] = NULL;
    ls->B2_blocks[0] = NULL;
    ls->C_blocks[0] = NULL;
}

/**
 * Find which block a state belongs to based on its non-zero pattern
 */
int btf_find_block(LS* ls, int* comb)
{
    int level = get_block_level(comb, ls->r);
    
    // Linear search through blocks to find matching pattern
    // In a production implementation, this could use a hash table
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (ls->block_levels[i] == level) {
            // Check if pattern matches
            int match = 1;
            for (int j = 0; j < ls->r; j++) {
                if ((comb[j] > 0) != (ls->block_patterns[i][j] > 0)) {
                    match = 0;
                    break;
                }
            }
            if (match) return i;
        }
    }
    
    // If no block found, create pattern for first block at this level
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (ls->block_levels[i] == level && ls->block_patterns[i][0] == 0) {
            // Initialize this block's pattern
            for (int j = 0; j < ls->r; j++) {
                ls->block_patterns[i][j] = (comb[j] > 0) ? 1 : 0;
            }
            return i;
        }
    }
    
    return 0; // Default to first block if no match
}

/**
 * Map local indices within a block to global indices
 */
int btf_local_to_global(LS* ls, int block_idx, int local_idx)
{
    return ls->block_positions[block_idx].start_row + local_idx;
}

/**
 * Free BTF block structures
 */
void btf_free_blocks(LS* ls)
{
    if (!ls) return;
    
    // Free X blocks
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (ls->X_blocks && ls->X_blocks[i]) {
            if (ls->X_blocks[i]->diag) {
                int size = ls->block_positions[i].rows;
                for (int j = 0; j < size; j++) {
                    if (ls->X_blocks[i]->diag[j]) {
                        for (int k = 0; k < size; k++) {
                            mpq_clear(ls->X_blocks[i]->diag[j][k]);
                        }
                        free(ls->X_blocks[i]->diag[j]);
                    }
                }
                free(ls->X_blocks[i]->diag);
            }
            if (ls->X_blocks[i]->lu_indices) {
                free(ls->X_blocks[i]->lu_indices);
            }
            free(ls->X_blocks[i]);
        }
    }
    
    // Free Y blocks
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (ls->Y_blocks && ls->Y_blocks[i]) {
            if (ls->Y_blocks[i]->nondiag) {
                if (ls->Y_blocks[i]->nondiag->coeff) {
                    for (int j = 0; j < ls->Y_blocks[i]->nondiag->nnz; j++) {
                        mpq_clear(ls->Y_blocks[i]->nondiag->coeff[j]);
                    }
                    free(ls->Y_blocks[i]->nondiag->coeff);
                }
                if (ls->Y_blocks[i]->nondiag->pos_row) free(ls->Y_blocks[i]->nondiag->pos_row);
                if (ls->Y_blocks[i]->nondiag->pos_col) free(ls->Y_blocks[i]->nondiag->pos_col);
                free(ls->Y_blocks[i]->nondiag);
            }
            free(ls->Y_blocks[i]);
        }
    }
    
    // Similar cleanup for B1, B2, C blocks
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (ls->B1_blocks && ls->B1_blocks[i]) {
            if (ls->B1_blocks[i]->nondiag) {
                if (ls->B1_blocks[i]->nondiag->coeff) {
                    for (int j = 0; j < ls->B1_blocks[i]->nondiag->nnz; j++) {
                        mpq_clear(ls->B1_blocks[i]->nondiag->coeff[j]);
                    }
                    free(ls->B1_blocks[i]->nondiag->coeff);
                }
                if (ls->B1_blocks[i]->nondiag->pos_row) free(ls->B1_blocks[i]->nondiag->pos_row);
                if (ls->B1_blocks[i]->nondiag->pos_col) free(ls->B1_blocks[i]->nondiag->pos_col);
                free(ls->B1_blocks[i]->nondiag);
            }
            free(ls->B1_blocks[i]);
        }
    }
    
    // Free block metadata
    if (ls->block_levels) free(ls->block_levels);
    if (ls->block_patterns) {
        for (int i = 0; i < ls->num_macro_blocks; i++) {
            if (ls->block_patterns[i]) free(ls->block_patterns[i]);
        }
        free(ls->block_patterns);
    }
    
    // Free block arrays
    if (ls->X_blocks) free(ls->X_blocks);
    if (ls->Y_blocks) free(ls->Y_blocks);
    if (ls->B1_blocks) free(ls->B1_blocks);
    if (ls->B2_blocks) free(ls->B2_blocks);
    if (ls->C_blocks) free(ls->C_blocks);
    if (ls->block_positions) free(ls->block_positions);
}