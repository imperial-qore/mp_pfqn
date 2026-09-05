#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <fpla.h>
#include "momf.h"

/* Floating-point counterpart of mom/perfindices.c.  Output formats are
 * kept byte-compatible with the exact solver so that the two can be
 * diffed directly; the exact num/den form (-e) has no meaning here and
 * falls back to the double form with a note on stderr. */
void perfindices(qnmodel* qnm, fp_vec_t G, fp_vec_t Gk, bool verbose_output, bool log_output, bool normconst_output, bool normconst_g_output, bool throughput_output, bool queue_output, bool debug_output, bool bounds_output, mpz_t scale_factor)
{
	int r;
	fp_t tmp; fp_init(tmp);
	fp_t tmp2; fp_init(tmp2);

	int Ntot = 0;
	for (r = 0; r < qnm->R; r++) Ntot += qnm->N[r];

	/* Undo the integer scaling introduced by the perturbation option. */
	fp_t G_scaled; fp_init(G_scaled);
	fp_set(G_scaled, G[0]);
	if (mpz_cmp_ui(scale_factor, 1) > 0)
	{
		mpz_t scale_power; mpz_init(scale_power);
		mpz_pow_ui(scale_power, scale_factor, Ntot);
		fp_t divisor; fp_init(divisor);
		fp_set_z(divisor, scale_power);
		fp_div(G_scaled, G_scaled, divisor);
		fp_clear(divisor);
		mpz_clear(scale_power);
	}

	/* log G is taken inside MPFR: at working precisions where G would
	 * overflow an IEEE double, mpfr_log still returns the right value,
	 * which is exactly the point of decoupling range from precision. */
	fp_t lg; fp_init(lg);
	mpfr_log(lg, G_scaled, MPFR_RNDN);
	double logG = fp_get_d(lg);
	fp_clear(lg);

	if (log_output)
	{
		printf("%.15e\n", logG);
		fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
		return;
	}
	else if (normconst_output)
	{
		fprintf(stderr, "note: -e has no exact form in the floating-point solver; printing G as a double\n");
		printf("%.15e\n", fp_get_d(G_scaled));
		fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
		return;
	}
	else if (normconst_g_output)
	{
		printf("%.15e\n", fp_get_d(G_scaled));
		fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
		return;
	}

	if (throughput_output)
	{
		for (r=1;r<=qnm->R;r++)
		{
			fp_div(G[r],G[r],G[0]);
			fp_t X_scaled; fp_init(X_scaled);
			fp_set(X_scaled, G[r]);
			if (mpz_cmp_ui(scale_factor, 1) > 0)
			{
				fp_t multiplier; fp_init(multiplier);
				fp_set_z(multiplier, scale_factor);
				fp_mul(X_scaled, X_scaled, multiplier);
				fp_clear(multiplier);
			}
			printf("%.15e\n", fp_get_d(X_scaled));
			fp_clear(X_scaled);
		}
		fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
		return;
	}

	if (queue_output)
	{
		for (int k=1; k<=qnm->M; k++)
		{
			for (r=1; r<=qnm->R; r++)
			{
				mpz_t original_L_value; mpz_init(original_L_value);
				mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
				if (mpz_cmp_ui(original_L_value, 0) == 0)
					fp_set_ui(tmp2, 0, 1);
				else
				{
					fp_set_z(tmp, qnm->L[k-1][r-1]);
					fp_mul(tmp2, Gk[(k-1)*(qnm->R+1)+r], tmp);
					fp_div(tmp2, tmp2, G[0]);
					fp_t mult; fp_init(mult);
					fp_set_si(mult, qnm->mi[k-1], 1);
					fp_mul(tmp2, tmp2, mult);
					fp_clear(mult);
				}
				mpz_clear(original_L_value);
				printf("%.15e", fp_get_d(tmp2));
				if (r < qnm->R) printf(" ");
			}
			printf("\n");
		}
		fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
		return;
	}

	/* default verbose output */
	if (bounds_output)
	{
		printf("========== Performance Metrics (Bounds) ==========\n");
		printf("Note: Bounds computation using dual perturbation is in development.\n");
		printf("Currently showing single perturbation results.\n");
	}
	else
		printf("========== Performance Metrics ==========\n");

	printf("G = %.15e\n", fp_get_d(G_scaled));
	printf("log(G) = %.15e\n", logG);

	printf("\n");
	printf("X (throughputs):\n");
	for (r=1;r<=qnm->R;r++)
	{
		fp_div(G[r],G[r],G[0]);
		fp_t X_scaled; fp_init(X_scaled);
		fp_set(X_scaled, G[r]);
		if (mpz_cmp_ui(scale_factor, 1) > 0)
		{
			fp_t multiplier; fp_init(multiplier);
			fp_set_z(multiplier, scale_factor);
			fp_mul(X_scaled, X_scaled, multiplier);
			fp_clear(multiplier);
		}
		printf("X[%d] = %.15e\n", r, fp_get_d(X_scaled));
		fp_clear(X_scaled);
	}

	printf("\nQ (mean queue lengths):\n");
	for (int k=1; k<=qnm->M; k++)
	{
		printf("Q[%d] =", k);
		fp_t total_q; fp_init(total_q); fp_set_ui(total_q, 0, 1);
		for (r=1; r<=qnm->R; r++)
		{
			mpz_t original_L_value; mpz_init(original_L_value);
			mpz_tdiv_q(original_L_value, qnm->L[k-1][r-1], scale_factor);
			if (mpz_cmp_ui(original_L_value, 0) == 0)
				fp_set_ui(tmp2, 0, 1);
			else
			{
				fp_set_z(tmp, qnm->L[k-1][r-1]);
				fp_mul(tmp2, Gk[(k-1)*(qnm->R+1)+r], tmp);
				fp_div(tmp2, tmp2, G[0]);
				fp_t mult; fp_init(mult);
				fp_set_si(mult, qnm->mi[k-1], 1);
				fp_mul(tmp2, tmp2, mult);
				fp_clear(mult);
			}
			mpz_clear(original_L_value);
			fp_add(total_q, total_q, tmp2);
			printf("\t%.15e", fp_get_d(tmp2));
		}
		printf("\t(total: %.15e)\n", fp_get_d(total_q));
		fp_clear(total_q);
	}

	if (bounds_output)
		printf("==================================================\n");
	else
		printf("=========================================\n");

	if (REFINE_ITERS > 0)
		printf("Refinement: %d sweeps over %d block solves (%.2f avg), working prec %ld bits, residual prec %ld bits\n",
		       REFINE_SWEEPS, REFINE_SOLVES,
		       REFINE_SOLVES ? (double)REFINE_SWEEPS/REFINE_SOLVES : 0.0,
		       (long)FP_WPREC, (long)FP_RPREC);
	else
		printf("Refinement: off, working prec %ld bits\n", (long)FP_WPREC);

	extern double t0, t1;
	t1 = CPUTIME;
	printf("Elapsed time (MoM-fp): %.6f s\n", t1-t0);

	fp_clear(tmp); fp_clear(tmp2); fp_clear(G_scaled);
}
