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

void printmodel_with_perturbation(qnmodel* qn, int perturbation_digit, long scale_factor, int perturbation_seed, mpz_t** original_L, mpz_t* original_Z)
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
	
	// Print think times (Z) with perturbation info if applicable
	printf("          Z[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++) {
		if (perturbation_digit > 0) {
			gmp_printf("%12Zd",original_Z[r-1]);
			// Phantom spacing to match L line formatting  
			printf("            ");  // 12 spaces to match " eps=X.0e-0X" format (12 chars)
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
				// Calculate and display perturbation with better randomization
				// Ensure all entries get perturbation, especially zeros
				long hash_val = perturbation_seed * 1103515245 + (m-1) * 12345 + (r-1) * 67891;
				long perturb = (abs(hash_val) % 9) + 1;  // 1-9 range, always positive
				double perturb_value = (double)perturb / (double)scale_factor;
				printf(" eps=%.1e", perturb_value);
			} else {
				gmp_printf("%12Zd",qn->L[m-1][r-1]);
			}
		}
		printf("\n");
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
