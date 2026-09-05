#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "util.h" 

mpz_t ncktable[MAXNCKTABLE][MAXNCKTABLE]; /* nchoosek table */
int ncktable_initialized = 0;

static void ensure_table_init() {
    if (!ncktable_initialized) {
        for (int i = 0; i < MAXNCKTABLE; i++) {
            for (int j = 0; j < MAXNCKTABLE; j++) {
                mpz_init(ncktable[i][j]);
            }
        }
        ncktable_initialized = 1;
    }
}

long long int nck(int N, int K)
{
    ensure_table_init();
    
    // Use symmetry: nck(n,k) = nck(n,n-k)
    // This avoids overflow for large k values
    if (K > N - K) {
        K = N - K;
    }
    
    if (N < 0 || K < 0 || K > N) {
        return 0;
    }
    
    if (K == 0 || K == N) {
        return 1;
    }
    
    if (N < MAXNCKTABLE && K < MAXNCKTABLE && mpz_cmp_ui(ncktable[N][K], 0) != 0) {
        // Check if the value fits in long long int
        if (mpz_fits_slong_p(ncktable[N][K])) {
            return mpz_get_si(ncktable[N][K]);
        } else {
            fprintf(stderr, "ERROR: nck(%d,%d) = ", N, K);
            mpz_out_str(stderr, 10, ncktable[N][K]);
            fprintf(stderr, " overflows long long int (max=%lld)\n", LLONG_MAX);
            fprintf(stderr, "This model has dimensions that are too large to handle.\n");
            exit(1);
        }
    }
    
    // Use exact arithmetic for computation
    mpz_t result, numerator, denominator;
    mpz_init(result);
    mpz_init(numerator);
    mpz_init(denominator);
    
    mpz_set_ui(numerator, 1);
    mpz_set_ui(denominator, 1);
    
    // Compute C(n,k) = n * (n-1) * ... * (n-k+1) / (k * (k-1) * ... * 1)
    for (int i = 0; i < K; i++) {
        mpz_mul_ui(numerator, numerator, N - i);
        mpz_mul_ui(denominator, denominator, i + 1);
    }
    
    mpz_divexact(result, numerator, denominator);
    
    // Cache the result
    if (N < MAXNCKTABLE && K < MAXNCKTABLE) {
        mpz_set(ncktable[N][K], result);
    }
    
    long long int ret;
    if (mpz_fits_slong_p(result)) {
        ret = mpz_get_si(result);
    } else {
        fprintf(stderr, "ERROR: nck(%d,%d) = ", N, K);
        mpz_out_str(stderr, 10, result);
        fprintf(stderr, " overflows long long int (max=%lld)\n", LLONG_MAX);
        fprintf(stderr, "This model has dimensions that are too large to handle.\n");
        exit(1);
    }
    
    mpz_clear(result);
    mpz_clear(numerator);
    mpz_clear(denominator);
    
    return ret;
}