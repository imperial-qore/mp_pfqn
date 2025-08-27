#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include "mom.h"
#ifdef MPFR
#include "mpfr.h"
#endif
#include <gmpla.h>

void perfindices(qnmodel* qnm, mpq_vec_t G, mpq_vec_t Gk, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor)
{
	int r;
	mpq_t tmp; mpq_init(tmp);
	mpq_t tmp2; mpq_init(tmp2);
	mpf_t fval; mpf_init(fval);
	
	// Compute total population for G scaling
	int Ntot = 0;
	for (r = 0; r < qnm->R; r++) {
		Ntot += qnm->N[r];
	}
	
	// Apply scaling correction to G if perturbation was used
	mpq_t G_scaled; mpq_init(G_scaled);
	mpq_set(G_scaled, G[0]);
	
	if (mpz_cmp_ui(scale_factor, 1) > 0) {
		// G needs to be divided by scale_factor^Ntot
		mpz_t scale_power;
		mpz_init(scale_power);
		mpz_pow_ui(scale_power, scale_factor, Ntot);
		mpq_t divisor;
		mpq_init(divisor);
		mpq_set_z(divisor, scale_power);
		mpq_div(G_scaled, G_scaled, divisor);
		mpq_clear(divisor);
		mpz_clear(scale_power);
	}
	
	mpf_set_q(fval, G_scaled);
	double logG = log(mpf_get_d(fval));
	
	if (log_output) {
		printf("%.15e\n", logG);
		mpq_clear(tmp);
		mpq_clear(tmp2);
		mpq_clear(G_scaled);
		mpf_clear(fval);
		return;
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
		mpq_clear(tmp);
		mpq_clear(tmp2);
		mpq_clear(G_scaled);
		mpf_clear(fval);
		return;
	} else if (normconst_g_output) {
		printf("%.15e\n", mpf_get_d(fval));
		mpq_clear(tmp);
		mpq_clear(tmp2);
		mpq_clear(G_scaled);
		mpf_clear(fval);
		return;
	}
	
	// Handle throughput-only output
	if (throughput_output) {
		for (r=1;r<=qnm->R;r++)
		{
			mpq_div(G[r],G[r],G[0]);
			
			// Apply scaling correction: X[r] needs to be multiplied by scale_factor
			mpq_t X_scaled; mpq_init(X_scaled);
			mpq_set(X_scaled, G[r]);
			
			if (mpz_cmp_ui(scale_factor, 1) > 0) {
				mpq_t multiplier;
				mpq_init(multiplier);
				mpq_set_z(multiplier, scale_factor);
				mpq_mul(X_scaled, X_scaled, multiplier);
				mpq_clear(multiplier);
			}
			
			mpf_set_q(fval, X_scaled);
			printf("%.15e\n", mpf_get_d(fval));
			mpq_clear(X_scaled);
		}
		mpq_clear(tmp);
		mpq_clear(tmp2);
		mpq_clear(G_scaled);
		mpf_clear(fval);
		return;
	}
	
	// Handle queue-length-only output
	if (queue_output) {
		for (int k=1; k<=qnm->M; k++)
		{
			for (r=1; r<=qnm->R; r++)
			{
				// Check if the original demand was 0.0
				mpz_t original_L_value;
				mpz_init(original_L_value);
				mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
				
				if (mpz_cmp_ui(original_L_value, 0) == 0) {
					// Original demand was 0, so Q should be 0
					mpq_set_ui(tmp2, 0, 1);
				} else {
					mpq_set_z(tmp, qnm->L[k-1][r-1]);
					mpq_mul(tmp2, Gk[(k-1)*(qnm->R+1)+r], tmp);
					mpf_set_q(fval, tmp2);
					mpq_div(tmp2, tmp2, G[0]);
					
					// Multiply by station multiplicity
					mpq_t mult; mpq_init(mult);
					mpq_set_si(mult, qnm->mi[k-1], 1);
					mpq_mul(tmp2, tmp2, mult);
					mpq_clear(mult);
				}
				
				mpz_clear(original_L_value);
				
				mpf_set_q(fval, tmp2);
				printf("%.15e", mpf_get_d(fval));
				if (r < qnm->R) {
					printf(" ");
				}
			}
			printf("\n");
		}
		mpq_clear(tmp);
		mpq_clear(tmp2);
		mpq_clear(G_scaled);
		mpf_clear(fval);
		return;
	}
	
	// Default verbose output
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		if (bounds_output) {
			printf("========== Performance Metrics (Bounds) ==========\n");
			printf("Note: Bounds computation using dual perturbation is in development.\n");
			printf("Currently showing single perturbation results.\n");
		} else {
			printf("========== Performance Metrics ==========\n");
		}
		if (debug_output) {
			// Print exact rational G, omit log(G)
			gmp_printf("G = %Qd\n", G_scaled);
		} else {
			// Print double G and log(G)
			printf("G = %.15e\n", mpf_get_d(fval));
			printf("log(G) = %.15e\n", logG);
		}
	}
	
	if (!log_output && !normconst_output && !normconst_g_output && !throughput_output && !queue_output) {
		printf("\n");
		printf("X (throughputs):\n");
		for (r=1;r<=qnm->R;r++)
		{
			mpq_div(G[r],G[r],G[0]);
			
			// Apply scaling correction: X[r] needs to be multiplied by scale_factor
			mpq_t X_scaled; mpq_init(X_scaled);
			mpq_set(X_scaled, G[r]);
			
			if (mpz_cmp_ui(scale_factor, 1) > 0) {
				mpq_t multiplier;
				mpq_init(multiplier);
				mpq_set_z(multiplier, scale_factor);
				mpq_mul(X_scaled, X_scaled, multiplier);
				mpq_clear(multiplier);
			}
			
			if (debug_output) {
				// Print exact rational throughputs
				gmp_printf("X[%d] = %Qd\n", r, X_scaled);
			} else {
				// Print double throughputs
				mpf_set_q(fval, X_scaled);
				printf("X[%d] = %.15e\n", r, mpf_get_d(fval));
			}
			mpq_clear(X_scaled);
		}
		
		printf("\nQ (mean queue lengths):\n");
		
		
		for (int k=1; k<=qnm->M; k++)
		{
			printf("Q[%d] =", k);
			mpq_t total_q; mpq_init(total_q); mpq_set_ui(total_q, 0, 1);
			for (r=1; r<=qnm->R; r++)
			{
				// Check if the original demand was 0.0
				mpz_t original_L_value;
				mpz_init(original_L_value);
				mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
				
				if (mpz_cmp_ui(original_L_value, 0) == 0) {
					// Original demand was 0, so Q should be 0
					mpq_set_ui(tmp2, 0, 1);
				} else {
					mpq_set_z(tmp, qnm->L[k-1][r-1]);
					mpq_mul(tmp2, Gk[(k-1)*(qnm->R+1)+r], tmp);
					mpq_div(tmp2, tmp2, G[0]);
					
					// Multiply by station multiplicity
					mpq_t mult; mpq_init(mult);
					mpq_set_si(mult, qnm->mi[k-1], 1);
					mpq_mul(tmp2, tmp2, mult);
					mpq_clear(mult);
				}
				
				mpz_clear(original_L_value);
				
				mpq_add(total_q, total_q, tmp2);
				
				if (debug_output) {
					// Print exact rational queue lengths
					printf("\t");
					gmp_printf("%Qd", tmp2);
				} else {
					// Print double queue lengths
					mpf_set_q(fval, tmp2);
					printf("\t%.15e", mpf_get_d(fval));
				}
			}
			
			if (debug_output) {
				printf("\t(total: ");
				gmp_printf("%Qd", total_q);
				printf(")\n");
			} else {
				mpf_set_q(fval, total_q);
				printf("\t(total: %.15e)\n", mpf_get_d(fval));
			}
			mpq_clear(total_q);
		}
		
		if (bounds_output) {
			printf("==================================================\n");
		} else {
			printf("=========================================\n");
		}
		extern double t0, t1;
		t1 = CPUTIME;
		printf("Elapsed time (MoM): %.6f s\n", t1-t0);
	}
	mpq_clear(tmp);
	mpq_clear(tmp2);
	mpq_clear(G_scaled);
	mpf_clear(fval);
}
