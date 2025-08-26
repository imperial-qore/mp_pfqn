#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"

void perfindices(qnmodel* qnm, mpq_vec_t g, mpq_t* G, mpq_t** Gk,
                 bool verbose_output, bool log_output, bool normconst_output, 
                 bool normconst_g_output, bool throughput_output, bool queue_output, 
                 long scale_factor) {
    
    mpf_t fval;
    mpf_init(fval);
    
    
    // G(N) is stored at position cardGk in the g vector
    long int finalCardGk = nck(qnm->M+qnm->R-1,qnm->M)*qnm->M;
    mpq_t G_total;
    mpq_init(G_total);
    mpq_set(G_total, g[finalCardGk]);
    
    // Apply scaling correction to normalizing constant
    mpq_t G_scaled;
    mpq_init(G_scaled);
    mpq_set(G_scaled, G_total);
    
    if (scale_factor > 1) {
        mpz_t divisor;
        mpz_init(divisor);
        mpz_ui_pow_ui(divisor, scale_factor, qnm->R);
        mpz_t num, den;
        mpz_init(num);
        mpz_init(den);
        mpq_get_num(num, G_scaled);
        mpq_get_den(den, G_scaled);
        mpz_mul(den, den, divisor);
        mpq_set_num(G_scaled, num);
        mpq_set_den(G_scaled, den);
        mpq_canonicalize(G_scaled);
        mpz_clear(divisor);
        mpz_clear(num);
        mpz_clear(den);
    }
    
    mpf_set_q(fval, G_scaled);
    double logG = log(mpf_get_d(fval));
    
    // Handle different output modes
    if (log_output) {
        printf("%.15e\n", logG);
    } else if (normconst_output) {
        // Print exact numerator and denominator
        mpz_t num, den;
        mpz_init(num);
        mpz_init(den);
        mpq_get_num(num, G_scaled);
        mpq_get_den(den, G_scaled);
        gmp_printf("%Zd\n", num);
        gmp_printf("%Zd\n", den);
        mpz_clear(num);
        mpz_clear(den);
    } else if (normconst_g_output) {
        // Print normalizing constant as double
        printf("%.15e\n", mpf_get_d(fval));
    } else if (throughput_output) {
        // Print throughputs, one per row
        mpq_t X_r;
        mpq_init(X_r);
        for (int r = 1; r <= qnm->R; r++) {
            mpq_div(X_r, G[r-1], G_total);
            
            // Apply scaling correction
            mpq_t X_scaled;
            mpq_init(X_scaled);
            mpq_set(X_scaled, X_r);
            
            if (scale_factor > 1) {
                mpq_t multiplier;
                mpq_init(multiplier);
                mpq_set_ui(multiplier, scale_factor, 1);
                mpq_mul(X_scaled, X_scaled, multiplier);
                mpq_clear(multiplier);
            }
            
            mpf_set_q(fval, X_scaled);
            printf("%.15e\n", mpf_get_d(fval));
            mpq_clear(X_scaled);
        }
        mpq_clear(X_r);
    } else if (queue_output) {
        // Print queue lengths, all classes for same queue on same row
        mpq_t Q_kr, tmp, tmp2;
        mpq_init(Q_kr);
        mpq_init(tmp);
        mpq_init(tmp2);
        for (int k = 1; k <= qnm->M; k++) {
            for (int r = 1; r <= qnm->R; r++) {
                // Check if the original demand was 0.0
                mpz_t original_L_value;
                mpz_init(original_L_value);
                mpz_tdiv_q_ui(original_L_value, qnm->L[k-1][r-1], scale_factor);
                
                if (mpz_cmp_ui(original_L_value, 0) == 0) {
                    // Original demand was 0, so Q should be 0
                    mpq_set_ui(Q_kr, 0, 1);
                } else {
                    mpq_set_z(tmp, qnm->L[k-1][r-1]);
                    mpq_mul(tmp2, Gk[k-1][r-1], tmp);
                    mpq_div(Q_kr, tmp2, G_total);
                }
                
                mpz_clear(original_L_value);
                
                mpf_set_q(fval, Q_kr);
                printf("%.15e", mpf_get_d(fval));
                if (r < qnm->R) printf(" ");
            }
            printf("\n");
        }
        mpq_clear(Q_kr);
        mpq_clear(tmp);
        mpq_clear(tmp2);
    } else {
        // Verbose output
        printf("========== Performance Metrics ==========\n");
        
        printf("G = %.15e\n", mpf_get_d(fval));
        printf("log(G) = %.15e\n", logG);
        
        // Compute exact throughputs using stored marginal normalizing constants
        printf("\nX (throughputs):\n");
        mpq_t X_r;
        mpq_init(X_r);
        
        for (int r = 1; r <= qnm->R; r++) {
            mpq_div(X_r, G[r-1], G_total);
            
            // Apply scaling correction
            mpq_t X_scaled;
            mpq_init(X_scaled);
            mpq_set(X_scaled, X_r);
            
            if (scale_factor > 1) {
                mpq_t multiplier;
                mpq_init(multiplier);
                mpq_set_ui(multiplier, scale_factor, 1);
                mpq_mul(X_scaled, X_scaled, multiplier);
                mpq_clear(multiplier);
            }
            
            mpf_set_q(fval, X_scaled);
            printf("X[%d] = %.15e\n", r, mpf_get_d(fval));
            mpq_clear(X_scaled);
        }
        
        // Compute queue lengths
        printf("\nQ (mean queue lengths):\n");
        
        mpq_t Q_kr, tmp, tmp2, total_q;
        mpq_init(Q_kr);
        mpq_init(tmp);
        mpq_init(tmp2);
        
        for (int k = 1; k <= qnm->M; k++) {
            printf("Q[%d] =", k);
            mpq_init(total_q);
            mpq_set_ui(total_q, 0, 1);
            
            for (int r = 1; r <= qnm->R; r++) {
                // Check if the original demand was 0.0
                mpz_t original_L_value;
                mpz_init(original_L_value);
                mpz_tdiv_q_ui(original_L_value, qnm->L[k-1][r-1], scale_factor);
                
                if (mpz_cmp_ui(original_L_value, 0) == 0) {
                    // Original demand was 0, so Q should be 0
                    mpq_set_ui(Q_kr, 0, 1);
                } else {
                    // Q[k,r] = L[k,r] * Gk^k(N-e_r) / G(N)
                    mpq_set_z(tmp, qnm->L[k-1][r-1]);
                    mpq_mul(tmp2, Gk[k-1][r-1], tmp);
                    mpq_div(Q_kr, tmp2, G_total);
                    
                    // Multiply by the multiplicity mi[k-1] for station k
                    if (qnm->mi[k-1] > 1) {
                        mpq_t mi_q;
                        mpq_init(mi_q);
                        mpq_set_ui(mi_q, qnm->mi[k-1], 1);
                        mpq_mul(Q_kr, Q_kr, mi_q);
                        mpq_clear(mi_q);
                    }
                }
                
                mpz_clear(original_L_value);
                
                mpq_add(total_q, total_q, Q_kr);
                
                mpf_set_q(fval, Q_kr);
                printf("\t%.15e", mpf_get_d(fval));
            }
            
            mpf_set_q(fval, total_q);
            printf("\t(total: %.15e)\n", mpf_get_d(fval));
            mpq_clear(total_q);
        }
        
        mpq_clear(Q_kr);
        mpq_clear(tmp);
        mpq_clear(tmp2);
        
        printf("=========================================\n");
        extern double t0, t1;
        t1 = CPUTIME;
        printf("Elapsed time (CoMoM): %g s\n", t1-t0);
        
        mpq_clear(X_r);
    }
    
    mpq_clear(G_total);
    mpq_clear(G_scaled);
    mpf_clear(fval);
}