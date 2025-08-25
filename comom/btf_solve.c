#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"

/**
 * BTF solver using backward substitution through block hierarchy
 */
mpq_vec_t btf_solve(LS* ls, mpq_vec_t b, int n)
{
    if (!ls) {
        return NULL;
    }
    
    int M = ls->m;
    int R = ls->r;
    
    if (M <= 0 || R <= 0) {
        return NULL;
    }
    
    long int cardG = nck(M + R - 1, M);
    long int cardGk = cardG * M;
    
    // For compatibility, use traditional solving if BTF blocks not initialized
    if (!ls->X_blocks || ls->num_macro_blocks == 0) {
        // Fallback to traditional solving
        mpq_vec_t G = (mpq_vec_t)mpq_vec(cardG, 0, 1);
        mpq_vec_t b1 = (mpq_vec_t)mpq_vec(cardGk, 0, 1);
        mpq_vec_t b1b = (mpq_vec_t)mpq_vec(cardGk, 0, 1);
        
        // Compute G = B2 * g
        mpq_mspvecmul(G, ls->B2, b);
        
        // Compute b1 = B1 * g - A12 * G
        mpq_mspvecmul(b1, ls->B1, b);
        mpq_mspvecmul(b1b, ls->A12, G);
        for (int i = 0; i < cardGk; i++) {
            mpq_sub(b1[i], b1[i], b1b[i]);
        }
        
        // LU decomposition (do it each time since A11 changes)
        int* lu_indices = mpq_ludcmp(ls->A11, cardGk);
        
        // Solve A11 * Gk = b1
        if (mpq_lubksb(ls->A11, b1, cardGk, lu_indices) < 0) {
            // System is singular
            free(lu_indices);
            for (int i = 0; i < cardG; i++) mpq_clear(G[i]);
            for (int i = 0; i < cardGk; i++) mpq_clear(b1[i]);
            for (int i = 0; i < cardGk; i++) mpq_clear(b1b[i]);
            free(G);
            free(b1);
            free(b1b);
            return NULL;
        }
        free(lu_indices);
        
        // Create result vector
        mpq_vec_t result = (mpq_vec_t)mpq_vec(cardG + cardGk, 0, 1);
        
        // Copy Gk and G into result
        for (int i = 0; i < cardGk; i++) {
            mpq_set(result[i], b1[i]);
        }
        for (int i = 0; i < cardG; i++) {
            mpq_set(result[cardGk + i], G[i]);
        }
        
        // Clean up temporary vectors
        for (int i = 0; i < cardG; i++) mpq_clear(G[i]);
        for (int i = 0; i < cardGk; i++) mpq_clear(b1[i]);
        for (int i = 0; i < cardGk; i++) mpq_clear(b1b[i]);
        free(G);
        free(b1);
        free(b1b);
        
        return result;
    }
    
    // Full BTF solving with block decomposition
    mpq_vec_t solution = (mpq_vec_t)mpq_vec(cardG + cardGk, 0, 1);
    mpq_vec_t* block_solutions = (mpq_vec_t*)calloc(ls->num_macro_blocks, sizeof(mpq_vec_t));
    
    // Solve blocks in reverse order (bottom to top of BTF)
    for (int level = ls->H; level >= 1; level--) {
        // Process all blocks at this level
        for (int block_idx = 0; block_idx < ls->num_macro_blocks; block_idx++) {
            if (ls->block_levels[block_idx] != level) continue;
            
            int block_size = ls->block_positions[block_idx].rows;
            mpq_vec_t block_rhs = (mpq_vec_t)mpq_vec(block_size, 0, 1);
            
            // Extract relevant part of RHS for this block
            int start = ls->block_positions[block_idx].start_row;
            for (int i = 0; i < block_size && start + i < cardGk; i++) {
                mpq_set(block_rhs[i], b[start + i]);
            }
            
            // Adjust RHS based on solutions from lower levels
            if (level < ls->H && ls->Y_blocks[block_idx]) {
                // Subtract contributions from already solved blocks
                mpq_vec_t temp = (mpq_vec_t)mpq_vec(block_size, 0, 1);
                for (int lower_block = block_idx + 1; lower_block < ls->num_macro_blocks; lower_block++) {
                    if (block_solutions[lower_block]) {
                        // Apply Y_block coupling
                        // This would use sparse matrix-vector multiplication
                        // For now, simplified
                    }
                }
                for (int i = 0; i < block_size; i++) {
                    mpq_sub(block_rhs[i], block_rhs[i], temp[i]);
                    mpq_clear(temp[i]);
                }
                free(temp);
            }
            
            // Solve diagonal block system
            if (!ls->X_blocks[block_idx]->lu_indices) {
                ls->X_blocks[block_idx]->lu_indices = mpq_ludcmp(ls->X_blocks[block_idx]->diag, block_size);
            }
            
            block_solutions[block_idx] = (mpq_vec_t)mpq_vec(block_size, 0, 1);
            for (int i = 0; i < block_size; i++) {
                mpq_set(block_solutions[block_idx][i], block_rhs[i]);
            }
            
            if (mpq_lubksb(ls->X_blocks[block_idx]->diag, block_solutions[block_idx], 
                          block_size, ls->X_blocks[block_idx]->lu_indices) < 0) {
                // Block is singular
                for (int i = 0; i < block_size; i++) {
                    mpq_clear(block_rhs[i]);
                }
                free(block_rhs);
                // Clean up and return NULL
                for (int i = 0; i < ls->num_macro_blocks; i++) {
                    if (block_solutions[i]) {
                        int size = ls->block_positions[i].rows;
                        for (int j = 0; j < size; j++) {
                            mpq_clear(block_solutions[i][j]);
                        }
                        free(block_solutions[i]);
                    }
                }
                free(block_solutions);
                for (int i = 0; i < cardG + cardGk; i++) {
                    mpq_clear(solution[i]);
                }
                free(solution);
                return NULL;
            }
            
            // Copy block solution to global solution vector
            for (int i = 0; i < block_size && start + i < cardGk; i++) {
                mpq_set(solution[start + i], block_solutions[block_idx][i]);
            }
            
            // Clean up block RHS
            for (int i = 0; i < block_size; i++) {
                mpq_clear(block_rhs[i]);
            }
            free(block_rhs);
        }
    }
    
    // Compute G part using B2
    mpq_vec_t G = (mpq_vec_t)mpq_vec(cardG, 0, 1);
    mpq_mspvecmul(G, ls->B2, b);
    for (int i = 0; i < cardG; i++) {
        mpq_set(solution[cardGk + i], G[i]);
        mpq_clear(G[i]);
    }
    free(G);
    
    // Clean up block solutions
    for (int i = 0; i < ls->num_macro_blocks; i++) {
        if (block_solutions[i]) {
            int size = ls->block_positions[i].rows;
            for (int j = 0; j < size; j++) {
                mpq_clear(block_solutions[i][j]);
            }
            free(block_solutions[i]);
        }
    }
    free(block_solutions);
    
    return solution;
}