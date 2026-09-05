#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mvaldmx.h"

/*
 * Compute effective capacity terms for MVA-LD-MX.
 * Port of pfqn_mvaldmx_ec.m from LINE-dev.
 *
 * All arithmetic is exact rational (mpq_t).
 *
 * mu is [M][mu_cols] where mu_cols >= Nt+1 (the caller extends it).
 * C[k] = 1/mu[ist][k], using 0-based indexing.
 *
 * MATLAB uses 1-based mu indexing: mu(ist, k) for k=1..Nt+extensions.
 * In C we use 0-based: mu[ist][k] for k=0..mu_cols-1.
 * MATLAB C(ist, k) = 1/mu(ist, k) maps to C[k] = 1/mu[ist][k-1] in 0-based.
 * But since mu is already 0-based in C, C[k] = 1/mu[ist][k].
 *
 * MATLAB b(ist) is 1-based: first k where mu(ist,k)==mu(ist,end).
 * In C, b is 0-based: first k where mu[ist][k]==mu[ist][Nt-1].
 *
 * MATLAB E(ist, 1+n) for n=0..Nt -> C: E[ist][n] for n=0..Nt
 * MATLAB EC(ist, n) for n=1..Nt -> C: EC[ist][n-1] for n=1..Nt
 */

void mvaldmx_ec(mpq_t **EC, mpq_t **E, mpq_t **Eprime,
                mpq_t *lambda, mpz_t **L, mpq_t **mu,
                int M, int R, int Nt, int mu_cols)
{
	int ist, n, n0, j, r;
	mpq_t Lo, tmp, tmp2;
	mpq_init(Lo);
	mpq_init(tmp);
	mpq_init(tmp2);

	/* Allocate C matrix: C[ist][k] = 1/mu[ist][k] for k=0..mu_cols-1 */
	mpq_t **C = (mpq_t **)malloc(M * sizeof(mpq_t *));
	for (ist = 0; ist < M; ist++) {
		C[ist] = (mpq_t *)malloc(mu_cols * sizeof(mpq_t));
		for (int k = 0; k < mu_cols; k++) {
			mpq_init(C[ist][k]);
			if (mpq_sgn(mu[ist][k]) != 0) {
				mpq_set_ui(C[ist][k], 1, 1);
				mpq_div(C[ist][k], C[ist][k], mu[ist][k]);
			}
			/* else C[ist][k] = 0 (already initialized) */
		}
	}

	/* Compute b[ist]: first 0-based index k where mu[ist][k] == mu[ist][Nt-1].
	 * MATLAB: b(ist) = find(mu(ist,:)==mu(ist,end), 1) -- 1-based
	 * Note: "end" in MATLAB is the original Nt (before extension), i.e., mu(ist, Nt)
	 * which is mu[ist][Nt-1] in 0-based. But after extension mu_cols > Nt, and
	 * mu(ist,end) in the extended array is still the same value. The find uses
	 * the original mu before extension (size Nt). So b is found in the first Nt entries.
	 *
	 * Actually looking at the MATLAB more carefully:
	 *   Nt = size(mu,2) is computed AFTER mu(:,end+1:end+1+max(b)) extension.
	 *   Wait, no: Nt = size(mu,2) is done BEFORE the extension in pfqn_mvaldmx_ec.
	 *   mu is passed in with the extension already done by the caller.
	 *   Actually looking at pfqn_mvaldmx.m: mu(:,end+1) = mu(:,end) first,
	 *   then calls pfqn_mvaldmx_ec(lambda,D,mu).
	 *   In pfqn_mvaldmx_ec: Nt = size(mu,2) which is original Nt + 1.
	 *   Then mu is extended further: mu(:,end+1:end+1+max(b)).
	 *   b(ist) = find(mu(ist,:)==mu(ist,end), 1) -- this uses the Nt+1 size mu.
	 *
	 * In our C code: mu_cols is the number of columns passed in.
	 * The caller (mvaldmx_solve) will extend mu by 1 column before calling us.
	 * So mu_cols = Nt + 1 when we receive it, and Nt here is the original Nt.
	 *
	 * We'll set Nt_local = mu_cols (the number of entries per station).
	 * b[ist] is 1-based in MATLAB. We'll use 1-based here too to match MATLAB indexing.
	 *
	 * Then we extend C by max(b) more elements (matching the MATLAB extension).
	 */

	/* b[ist] = first 1-based index k where mu[ist][k-1] == mu[ist][mu_cols-1] */
	int *b = (int *)calloc(M, sizeof(int));
	for (ist = 0; ist < M; ist++) {
		b[ist] = 1; /* default */
		for (int k = 0; k < mu_cols; k++) {
			if (mpq_equal(mu[ist][k], mu[ist][mu_cols - 1])) {
				b[ist] = k + 1; /* 1-based */
				break;
			}
		}
	}

	int maxb = 0;
	for (ist = 0; ist < M; ist++)
		if (b[ist] > maxb) maxb = b[ist];

	/* In MATLAB: Nt = size(mu,2) after the first extension (Nt+1 columns).
	 * Then mu is extended further. In our code, mu_cols is already Nt+1.
	 * We set Nt_ec = mu_cols as the "Nt" used in this function. */
	int Nt_ec = mu_cols;

	/* Extend C to have enough columns: we need C[ist][k] for k up to
	 * at least Nt_ec-1+maxb or so. Let's compute the max index needed.
	 * In MATLAB: C(ist, n+n0) where n can be up to b-1 and n0 up to b-2.
	 * Max index = (b-1) + (b-2) = 2*b-3 (1-based), so 0-based = 2*b-4.
	 * Also need C(ist, n+n0+1) for F2prime, max = b-1 + b-2 + 1 = 2*b-2 (1-based).
	 * We also need up to Nt_ec (1-based) for EC computation.
	 * Safe upper bound: Nt_ec + maxb + 2. */
	int C_cols = Nt_ec + maxb + 2;
	/* Reallocate C with extended size */
	for (ist = 0; ist < M; ist++) {
		C[ist] = (mpq_t *)realloc(C[ist], C_cols * sizeof(mpq_t));
		for (int k = mu_cols; k < C_cols; k++) {
			mpq_init(C[ist][k]);
			/* Extended entries use mu[ist][mu_cols-1] (the saturation rate) */
			if (mpq_sgn(mu[ist][mu_cols - 1]) != 0) {
				mpq_set_ui(C[ist][k], 1, 1);
				mpq_div(C[ist][k], C[ist][k], mu[ist][mu_cols - 1]);
			}
		}
	}

	/* Now compute E, Eprime, EC per station.
	 * We use 1-based MATLAB indexing internally and map to 0-based arrays.
	 * MATLAB: E(ist, 1+n) for n=0..Nt, C(ist, k) for k=1..
	 * C array: C[ist][k-1] for 1-based k, i.e. C[ist][0] = 1/mu[ist][0] = C(ist,1).
	 *
	 * Actually, let's be careful. In MATLAB:
	 *   C = 1./mu  -- so C(ist, k) = 1/mu(ist, k), k=1..size(mu,2)
	 *   C(ist, b(ist)) means 1/mu(ist, b(ist))
	 *
	 * In C, our C array is 0-based: C[ist][k] = 1/mu[ist][k].
	 * MATLAB C(ist, k) = our C[ist][k-1].
	 *
	 * b(ist) in MATLAB is 1-based. Let's define b_m = b[ist] (1-based).
	 * C_b = C[ist][b_m - 1] = MATLAB C(ist, b(ist)).
	 */

	/* Allocate local arrays for E1, E2, E3, F2, F3, F2prime, E2prime per station */
	for (ist = 0; ist < M; ist++) {
		int b_m = b[ist]; /* 1-based saturation level */
		/* C_b = C[ist][b_m - 1] */
		/* Lo(ist) = lambda * D(ist,:)' = sum_r lambda[r] * L[ist][r] */
		mpq_set_ui(Lo, 0, 1);
		for (r = 0; r < R; r++) {
			if (mpq_sgn(lambda[r]) != 0 && mpz_sgn(L[ist][r]) != 0) {
				mpq_set_z(tmp, L[ist][r]);
				mpq_mul(tmp, tmp, lambda[r]);
				mpq_add(Lo, Lo, tmp);
			}
		}

		/* Allocate E1, E2, E3, E2prime arrays [0..Nt_ec] (Nt_ec+1 entries) */
		mpq_t *E1 = (mpq_t *)malloc((Nt_ec + 1) * sizeof(mpq_t));
		mpq_t *E2 = (mpq_t *)malloc((Nt_ec + 1) * sizeof(mpq_t));
		mpq_t *E3 = (mpq_t *)malloc((Nt_ec + 1) * sizeof(mpq_t));
		mpq_t *E2prime = (mpq_t *)malloc((Nt_ec + 1) * sizeof(mpq_t));
		for (n = 0; n <= Nt_ec; n++) {
			mpq_init(E1[n]);
			mpq_init(E2[n]);
			mpq_init(E3[n]);
			mpq_init(E2prime[n]);
		}

		/* F2[n][n0] for n=0..Nt_ec, n0=0..b_m-2 */
		/* F3[n][n0] for n=0..Nt_ec, n0=0..b_m-2 */
		/* F2prime[n][n0] for n=0..Nt_ec, n0=0..b_m-2 */
		int f_cols = (b_m >= 2) ? (b_m - 1) : 0; /* number of n0 values: 0..b_m-2 */
		mpq_t **F2 = NULL, **F3 = NULL, **F2prime_arr = NULL;
		if (f_cols > 0) {
			F2 = (mpq_t **)malloc((Nt_ec + 1) * sizeof(mpq_t *));
			F3 = (mpq_t **)malloc((Nt_ec + 1) * sizeof(mpq_t *));
			F2prime_arr = (mpq_t **)malloc((Nt_ec + 1) * sizeof(mpq_t *));
			for (n = 0; n <= Nt_ec; n++) {
				F2[n] = (mpq_t *)malloc(f_cols * sizeof(mpq_t));
				F3[n] = (mpq_t *)malloc(f_cols * sizeof(mpq_t));
				F2prime_arr[n] = (mpq_t *)malloc(f_cols * sizeof(mpq_t));
				for (n0 = 0; n0 < f_cols; n0++) {
					mpq_init(F2[n][n0]);
					mpq_init(F3[n][n0]);
					mpq_init(F2prime_arr[n][n0]);
				}
			}
		}

		for (n = 0; n <= Nt_ec; n++) {
			/* MATLAB: n ranges 0..Nt, and checks n >= b(ist) vs n <= b(ist)-1.
			 * b(ist) is 1-based. So n >= b_m means the saturation regime.
			 * n < b_m means the non-saturated regime. */
			if (n >= b_m) {
				/* E(ist,1+n) = 1/(1 - Lo*C_b)^(n+1)
				 * Eprime(ist,1+n) = C_b * E(ist,1+n) */
				/* denom = 1 - Lo * C_b */
				mpq_t denom, base;
				mpq_init(denom);
				mpq_init(base);
				mpq_mul(denom, Lo, C[ist][b_m - 1]);
				mpq_set_ui(base, 1, 1);
				mpq_sub(denom, base, denom); /* denom = 1 - Lo*C_b */
				/* E = 1 / denom^(n+1) */
				/* Compute denom^(n+1) by repeated multiplication */
				mpq_t denom_pow;
				mpq_init(denom_pow);
				mpq_set_ui(denom_pow, 1, 1);
				for (int p = 0; p < n + 1; p++)
					mpq_mul(denom_pow, denom_pow, denom);
				mpq_set_ui(E[ist][n], 1, 1);
				mpq_div(E[ist][n], E[ist][n], denom_pow);

				mpq_mul(Eprime[ist][n], C[ist][b_m - 1], E[ist][n]);

				mpq_clear(denom);
				mpq_clear(base);
				mpq_clear(denom_pow);
			} else {
				/* n < b_m: compute E1, F2, E2, F3, E3, F2prime, E2prime */

				/* --- E1 --- */
				if (n == 0) {
					/* E1(1+0) = 1/(1-Lo*C_b) * prod_{j=1}^{b-1} C(j)/C_b
					 * MATLAB j=1..b-1, C(ist,j) = our C[ist][j-1]
					 * C_b = C[ist][b_m-1] */
					mpq_t denom_e1;
					mpq_init(denom_e1);
					mpq_mul(denom_e1, Lo, C[ist][b_m - 1]);
					mpq_set_ui(tmp, 1, 1);
					mpq_sub(denom_e1, tmp, denom_e1);
					mpq_set_ui(E1[n], 1, 1);
					mpq_div(E1[n], E1[n], denom_e1);
					for (j = 1; j <= b_m - 1; j++) {
						/* E1 *= C[ist][j-1] / C[ist][b_m-1] */
						mpq_div(tmp, C[ist][j - 1], C[ist][b_m - 1]);
						mpq_mul(E1[n], E1[n], tmp);
					}
					mpq_clear(denom_e1);
				} else {
					/* n > 0: E1(1+n) = 1/(1-Lo*C_b) * C_b/C(ist,n) * E1(1+(n-1))
					 * C(ist,n) = C[ist][n-1] */
					mpq_t denom_e1;
					mpq_init(denom_e1);
					mpq_mul(denom_e1, Lo, C[ist][b_m - 1]);
					mpq_set_ui(tmp, 1, 1);
					mpq_sub(denom_e1, tmp, denom_e1);

					mpq_set_ui(tmp, 1, 1);
					mpq_div(tmp, tmp, denom_e1);
					/* C_b / C(ist,n) = C[ist][b_m-1] / C[ist][n-1] */
					mpq_div(tmp2, C[ist][b_m - 1], C[ist][n - 1]);
					mpq_mul(tmp, tmp, tmp2);
					mpq_mul(E1[n], tmp, E1[n - 1]);
					mpq_clear(denom_e1);
				}

				/* --- F2, E2 --- */
				if (f_cols > 0) {
					for (n0 = 0; n0 <= b_m - 2; n0++) {
						if (n0 == 0) {
							mpq_set_ui(F2[n][n0], 1, 1);
						} else {
							/* F2(1+n, 1+n0) = (n+n0)/n0 * Lo * C(ist,n+n0) * F2(1+n, 1+(n0-1))
							 * C(ist, n+n0) = C[ist][n+n0-1] */
							mpq_set_ui(tmp, n + n0, n0);
							mpq_canonicalize(tmp);
							mpq_mul(tmp, tmp, Lo);
							mpq_mul(tmp, tmp, C[ist][n + n0 - 1]);
							mpq_mul(F2[n][n0], tmp, F2[n][n0 - 1]);
						}
					}
					/* E2(1+n) = sum(F2(1+n, 1+(0:b-2))) */
					mpq_set_ui(E2[n], 0, 1);
					for (n0 = 0; n0 <= b_m - 2; n0++)
						mpq_add(E2[n], E2[n], F2[n][n0]);
				} else {
					mpq_set_ui(E2[n], 0, 1);
				}

				/* --- F3, E3 --- */
				if (f_cols > 0) {
					for (n0 = 0; n0 <= b_m - 2; n0++) {
						if (n == 0 && n0 == 0) {
							/* F3(1,1) = prod_{j=1}^{b-1} C(ist,j)/C_b
							 * = prod_{j=1}^{b-1} C[ist][j-1] / C[ist][b_m-1] */
							mpq_set_ui(F3[n][n0], 1, 1);
							for (j = 1; j <= b_m - 1; j++) {
								mpq_div(tmp, C[ist][j - 1], C[ist][b_m - 1]);
								mpq_mul(F3[n][n0], F3[n][n0], tmp);
							}
						} else if (n > 0 && n0 == 0) {
							/* F3(1+n, 1+0) = C_b / C(ist,n) * F3(1+(n-1), 1+0)
							 * = C[ist][b_m-1] / C[ist][n-1] * F3[n-1][0] */
							mpq_div(tmp, C[ist][b_m - 1], C[ist][n - 1]);
							mpq_mul(F3[n][n0], tmp, F3[n - 1][0]);
						} else {
							/* F3(1+n, 1+n0) = (n+n0)/n0 * Lo * C_b * F3(1+n, 1+(n0-1))
							 * C_b = C[ist][b_m-1] */
							mpq_set_ui(tmp, n + n0, n0);
							mpq_canonicalize(tmp);
							mpq_mul(tmp, tmp, Lo);
							mpq_mul(tmp, tmp, C[ist][b_m - 1]);
							mpq_mul(F3[n][n0], tmp, F3[n][n0 - 1]);
						}
					}
					/* E3(1+n) = sum(F3(1+n, 1+(0:b-2))) */
					mpq_set_ui(E3[n], 0, 1);
					for (n0 = 0; n0 <= b_m - 2; n0++)
						mpq_add(E3[n], E3[n], F3[n][n0]);
				} else {
					mpq_set_ui(E3[n], 0, 1);
				}

				/* --- F2prime, E2prime --- */
				if (f_cols > 0) {
					for (n0 = 0; n0 <= b_m - 2; n0++) {
						if (n0 == 0) {
							/* F2prime(1+n, 1+0) = C(ist, n+1) = C[ist][n] */
							mpq_set(F2prime_arr[n][n0], C[ist][n]);
						} else {
							/* F2prime(1+n, 1+n0) = (n+n0)/n0 * Lo * C(ist, n+n0+1) * F2prime(1+n, 1+(n0-1))
							 * C(ist, n+n0+1) = C[ist][n+n0] */
							mpq_set_ui(tmp, n + n0, n0);
							mpq_canonicalize(tmp);
							mpq_mul(tmp, tmp, Lo);
							mpq_mul(tmp, tmp, C[ist][n + n0]);
							mpq_mul(F2prime_arr[n][n0], tmp, F2prime_arr[n][n0 - 1]);
						}
					}
					/* E2prime(1+n) = sum(F2prime(1+n, 1+(0:b-2))) */
					mpq_set_ui(E2prime[n], 0, 1);
					for (n0 = 0; n0 <= b_m - 2; n0++)
						mpq_add(E2prime[n], E2prime[n], F2prime_arr[n][n0]);
				} else {
					mpq_set_ui(E2prime[n], 0, 1);
				}

				/* E(ist, 1+n) = E1 + E2 - E3 */
				mpq_add(E[ist][n], E1[n], E2[n]);
				mpq_sub(E[ist][n], E[ist][n], E3[n]);

				/* Eprime */
				if (n < b_m - 1) {
					/* Eprime(ist,1+n) = C_b * E1(1+n) + E2prime(1+n) - C_b * E3(1+n) */
					mpq_mul(tmp, C[ist][b_m - 1], E1[n]);
					mpq_add(Eprime[ist][n], tmp, E2prime[n]);
					mpq_mul(tmp, C[ist][b_m - 1], E3[n]);
					mpq_sub(Eprime[ist][n], Eprime[ist][n], tmp);
				} else {
					/* n >= b_m - 1: Eprime(ist,1+n) = C_b * E(ist,1+n) */
					mpq_mul(Eprime[ist][n], C[ist][b_m - 1], E[ist][n]);
				}
			}
		}

		/* EC(ist, n) for n=1..Nt_ec: EC[ist][n-1] = C(ist,n) * E(ist,1+n) / E(ist,1+(n-1))
		 * C(ist,n) = C[ist][n-1] */
		for (n = 1; n <= Nt_ec; n++) {
			if (mpq_sgn(E[ist][n - 1]) != 0) {
				mpq_mul(EC[ist][n - 1], C[ist][n - 1], E[ist][n]);
				mpq_div(EC[ist][n - 1], EC[ist][n - 1], E[ist][n - 1]);
			}
		}

		/* Free local arrays */
		for (n = 0; n <= Nt_ec; n++) {
			mpq_clear(E1[n]);
			mpq_clear(E2[n]);
			mpq_clear(E3[n]);
			mpq_clear(E2prime[n]);
		}
		free(E1);
		free(E2);
		free(E3);
		free(E2prime);

		if (f_cols > 0) {
			for (n = 0; n <= Nt_ec; n++) {
				for (n0 = 0; n0 < f_cols; n0++) {
					mpq_clear(F2[n][n0]);
					mpq_clear(F3[n][n0]);
					mpq_clear(F2prime_arr[n][n0]);
				}
				free(F2[n]);
				free(F3[n]);
				free(F2prime_arr[n]);
			}
			free(F2);
			free(F3);
			free(F2prime_arr);
		}
	}

	/* Free C matrix */
	for (ist = 0; ist < M; ist++) {
		for (int k = 0; k < C_cols; k++)
			mpq_clear(C[ist][k]);
		free(C[ist]);
	}
	free(C);
	free(b);

	mpq_clear(Lo);
	mpq_clear(tmp);
	mpq_clear(tmp2);
}
