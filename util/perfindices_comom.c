#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

void perfindices_comom(qnmodel* qnm, mpq_t* g, combsrep* Dn)
{
    int r, s;
    mpq_t G, tmp, tmp2;
    mpq_t *Gs, *Gks, *X;
    mpf_t fval;
    
    mpq_init(G);
    mpq_init(tmp);
    mpq_init(tmp2);
    mpf_init(fval);
    
    // Allocate arrays for G_s, G^k_s, and X
    Gs = (mpq_t*)malloc((qnm->R+1)*sizeof(mpq_t));
    Gks = (mpq_t*)malloc((qnm->R+1)*sizeof(mpq_t));
    X = (mpq_t*)malloc((qnm->R+1)*sizeof(mpq_t));
    
    for(r=0; r<=qnm->R; r++) {
        mpq_init(Gs[r]);
        mpq_init(Gks[r]);
        mpq_init(X[r]);
    }
    
    // Get G (total normalizing constant) - this is g[N]
    int idx_N = hash(Dn, qnm->N, qnm->R);
    mpq_set(G, g[idx_N]);
    
    // Get G_s values - these are g[N-e_s]
    int *n_minus_es = (int*)malloc(qnm->R * sizeof(int));
    for(s=0; s<qnm->R; s++) {
        // Create N - e_s
        for(r=0; r<qnm->R; r++) {
            n_minus_es[r] = qnm->N[r];
        }
        n_minus_es[s]--;
        
        int idx = hash(Dn, n_minus_es, qnm->R);
        mpq_set(Gs[s+1], g[idx]);
    }
    
    // Get G^k_s values - these are g[N-2*e_s]
    for(s=0; s<qnm->R; s++) {
        // Create N - 2*e_s
        for(r=0; r<qnm->R; r++) {
            n_minus_es[r] = qnm->N[r];
        }
        n_minus_es[s] -= 2;
        
        if(n_minus_es[s] >= 0) {
            int idx = hash(Dn, n_minus_es, qnm->R);
            mpq_set(Gks[s+1], g[idx]);
        } else {
            mpq_set_si(Gks[s+1], 0, 1);
        }
    }
    
    // Calculate throughputs X_s = G_s / G
    for(s=1; s<=qnm->R; s++) {
        mpq_div(X[s], Gs[s], G);
    }
    
    // Print results
    printf("\n========== Performance Metrics ==========\n");
    
    // Print G with rational and double
    printf("G = ");
    gmp_printf("%Qd", G);
    mpf_set_q(fval, G);
    printf(" = %.15e\n", mpf_get_d(fval));
    
    // Print G_s values
    printf("G_s =");
    for(r=1; r<=qnm->R; r++) {
        printf("\t");
        gmp_printf("%Qd", Gs[r]);
    }
    printf("\n     ");
    for(r=1; r<=qnm->R; r++) {
        mpf_set_q(fval, Gs[r]);
        printf("\t%.15e", mpf_get_d(fval));
    }
    
    // Print G^k_s values
    printf("\nG^k_s =");
    for(r=1; r<=qnm->R; r++) {
        printf("\t");
        gmp_printf("%Qd", Gks[r]);
    }
    printf("\n       ");
    for(r=1; r<=qnm->R; r++) {
        mpf_set_q(fval, Gks[r]);
        printf("\t%.15e", mpf_get_d(fval));
    }
    
    // Print throughputs
    printf("\nX =\t");
    for(r=1; r<=qnm->R; r++) {
        gmp_printf("%Qd", X[r]);
        printf("\t");
    }
    printf("\n   \t");
    for(r=1; r<=qnm->R; r++) {
        mpf_set_q(fval, X[r]);
        printf("%.15e\t", mpf_get_d(fval));
    }
    
    // Calculate and print Q_11 = L_11 * G^k_1 / G
    if(qnm->M > 0 && qnm->R > 0) {
        mpq_set_si(tmp, qnm->L[0][0], 1);
        mpq_mul(tmp2, Gks[1], tmp);
        mpq_div(tmp2, tmp2, G);
        
        printf("\nQ_11 = ");
        gmp_printf("%Qd", tmp2);
        mpf_set_q(fval, tmp2);
        printf(" = %.15e", mpf_get_d(fval));
    }
    
    printf("\n=========================================\n");
    
    // Cleanup
    mpq_clear(G);
    mpq_clear(tmp);
    mpq_clear(tmp2);
    mpf_clear(fval);
    free(n_minus_es);
    
    for(r=0; r<=qnm->R; r++) {
        mpq_clear(Gs[r]);
        mpq_clear(Gks[r]);
        mpq_clear(X[r]);
    }
    free(Gs);
    free(Gks);
    free(X);
}