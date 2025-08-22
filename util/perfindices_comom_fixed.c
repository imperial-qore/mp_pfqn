#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

// Forward declaration for hash function from COMOM
extern int hash(combsrep* Dn, int *comb, int i);

void perfindices_comom(qnmodel* qnm, mpq_t* g, combsrep* Dn, bool full_output)
{
    int r;
    mpq_t G_total;
    mpf_t fval;
    
    mpq_init(G_total);
    mpf_init(fval);
    
    // Get G(N) - total normalizing constant
    // In COMOM, this is stored at position finalCardGk in the g vector
    long int finalCardGk = (long int)Dn->card * qnm->M;
    mpq_set(G_total, g[finalCardGk]);
    
    printf("\n========== Performance Metrics ==========\n");
    
    // Print G
    printf("G = ");
    if (full_output) {
        gmp_printf("%Qd", G_total);
        printf(" = ");
    }
    mpf_set_q(fval, G_total);
    printf("%.15e\n", mpf_get_d(fval));
    double logG = log(mpf_get_d(fval));
    printf("log(G) = %.15e\n", logG);
    
    // For throughputs, we need marginal normalizing constants G(N-e_r)
    // These are not directly available from the final g vector in the current COMOM implementation
    // The COMOM algorithm would need to be modified to store these values during computation
    
    printf("\nX (throughputs):\n");
    printf("Note: Throughput computation requires marginal normalizing constants.\n");
    printf("The current COMOM implementation does not store G(N-e_r) values.\n");
    printf("Full throughput computation requires algorithm modification.\n");
    
    printf("\nQ (queue lengths):\n");
    printf("Note: Queue length computation also requires marginal normalizing constants.\n");
    
    printf("=========================================\n");
    
    // Cleanup
    mpq_clear(G_total);
    mpf_clear(fval);
}