#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

void printmodel(qnmodel* qn)
{
	int m,r,mtot=0;
	for(m=1;m<=qn->M;m++)
		mtot+=qn->mi[m-1];
	printf("The queueing network model has %d queues (%d are replicas) and %d classes\n",mtot,mtot-qn->M,qn->R);
	printf("          N[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++)
		printf("%12d",qn->N[r-1]);
	printf("\n");
	printf("          Z[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++)
		gmp_printf("%12Zd",qn->Z[r-1]);
	printf("\n");
	for(m=1;m<=qn->M;m++)
	{
	 printf("mi=%d    L[%d,1:%d]:",qn->mi[m-1],m,qn->R);
	 for (r=1;r<=qn->R;r++)
	 {
		gmp_printf("%12Zd",qn->L[m-1][r-1]);
	 }
	 printf("\n");
	}
	printf("\n");
}

// Helper function to generate a permutation - same as in mom/main.c
static void generate_permutation_for_print(int* perm, int n, unsigned int seed) {
	// Initialize with 1 to n
	for(int i = 0; i < n; i++) {
		perm[i] = i + 1;
	}
	
	// Use a simple linear congruential generator for reproducible randomness
	unsigned int rand_state = seed;
	
	// Fisher-Yates shuffle
	for(int i = n - 1; i > 0; i--) {
		// Generate random number using LCG
		rand_state = rand_state * 1103515245 + 12345;
		int j = (rand_state / 65536) % (i + 1);
		
		// Swap elements i and j
		int temp = perm[i];
		perm[i] = perm[j];
		perm[j] = temp;
	}
}

void printmodel_with_perturbation(qnmodel* qn, int perturbation_digit, mpz_t scale_factor, int perturbation_seed, mpz_t** original_L, mpz_t* original_Z)
{
	int m,r,mtot=0;
	for(m=1;m<=qn->M;m++)
		mtot+=qn->mi[m-1];
	printf("The queueing network model has %d queues (%d are replicas) and %d classes\n",mtot,mtot-qn->M,qn->R);
	
	// Print populations (N) as integers - these are not perturbed
	printf("          N[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++) {
		printf("%12d",qn->N[r-1]);
		// Phantom spacing to match L line formatting
		if (perturbation_digit > 0) {
			printf("            ");  // 12 spaces to match " eps=X.0e-0X" format (12 chars)
		}
	}
	printf("\n");
	
	// Prepare permutation array
	int* perm = NULL;
	if (perturbation_digit > 0) {
		perm = (int*)calloc(qn->M + 1, sizeof(int));
	}
	
	// Print think times (Z) with perturbation info if applicable
	printf("          Z[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++) {
		if (perturbation_digit > 0) {
			gmp_printf("%12Zd",original_Z[r-1]);
			// Generate permutation for class r-1
			generate_permutation_for_print(perm, qn->M + 1, perturbation_seed + (r-1) * 1000);
			// Z gets the last element of the permutation
			double perturb_value = (double)perm[qn->M] / mpz_get_d(scale_factor);
			printf(" eps=%.1e", perturb_value);
		} else {
			gmp_printf("%12Zd",qn->Z[r-1]);
		}
	}
	printf("\n");
	
	// Print service demands (L) with perturbation info
	for(m=1;m<=qn->M;m++) {
		printf("mi=%d    L[%d,1:%d]:",qn->mi[m-1],m,qn->R);
		for (r=1;r<=qn->R;r++) {
			if (perturbation_digit > 0) {
				gmp_printf("%12Zd",original_L[m-1][r-1]);
				// Generate permutation for class r-1
				generate_permutation_for_print(perm, qn->M + 1, perturbation_seed + (r-1) * 1000);
				// L[m-1][r-1] gets perm[m-1]
				double perturb_value = (double)perm[m-1] / mpz_get_d(scale_factor);
				printf(" eps=%.1e", perturb_value);
			} else {
				gmp_printf("%12Zd",qn->L[m-1][r-1]);
			}
		}
		printf("\n");
	}
	
	if (perm != NULL) {
		free(perm);
	}
	printf("\n");
}

/*
int main(int argc, char** argv)
{
	printf("Opening %s\n",argv[1]);
	qnmodel* qn=(qnmodel*) readmodel(argv[1]);
	printmodel(qn);
	return 0;
}
*/
