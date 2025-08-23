#include <stdio.h>
#include "gmpla.h"
#include "util.h"
#include "assert.h"

qnmodel* readmodel(char* filename)
{
	qnmodel* qn=(qnmodel*)malloc(sizeof(qnmodel));
	int m;
	int r;

	/* open file*/
	FILE *f = fopen(filename,"r");

	if(!f)
	{
		printf("File %s cannot be opened\n", filename);
		exit(1);
	}	

	/* read R */
	if (fscanf(f,"%d\n",&qn->R) != 1) {
		printf("Error reading R from file\n");
		exit(1);
	}
	/* read N */
	qn->N=(int*)int_vec(qn->R,0);
	for (r=1;r<=qn->R;r++)
	{
		if (fscanf(f,"%d",&qn->N[r-1]) != 1) {
			printf("Error reading N[%d] from file\n", r-1);
			exit(1);
		}
	} 
	if (fscanf(f,"\n") != 0) {
		/* Newline consumed, continue */
	}
	/* read Z */
	qn->Z=(mpz_t*)malloc(qn->R * sizeof(mpz_t));
	for (r=1;r<=qn->R;r++)
	{
		mpz_init(qn->Z[r-1]);
		if (gmp_fscanf(f,"%Zd",&qn->Z[r-1]) != 1) {
			printf("Error reading Z[%d] from file\n", r-1);
			exit(1);
		}
	}
	if (fscanf(f,"\n") != 0) {
		/* Newline consumed, continue */
	}
	/* read M */
	if (fscanf(f,"%d\n",&qn->M) != 1) {
		printf("Error reading M from file\n");
		exit(1);
	}
	/* read mi and L */
	qn->L=(mpz_t**)malloc(qn->M * sizeof(mpz_t*));
	for(m=0; m<qn->M; m++) {
		qn->L[m] = (mpz_t*)malloc(qn->R * sizeof(mpz_t));
		for(r=0; r<qn->R; r++) {
			mpz_init(qn->L[m][r]);
		}
	}
	qn->mi=(int*)int_vec(qn->M,0);
	for(m=1;m<=qn->M;m++)
	{
 	 if (fscanf(f,"%d",&qn->mi[m-1]) != 1) {
		printf("Error reading mi[%d] from file\n", m-1);
		exit(1);
	 }
	 for (r=1;r<=qn->R;r++)
	 {
		if (gmp_fscanf(f,"%Zd",&qn->L[m-1][r-1]) != 1) {
			printf("Error reading L[%d][%d] from file\n", m-1, r-1);
			exit(1);
		}
	 }
	 if (fscanf(f,"\n") != 0) {
		/* Newline consumed, continue */
	 }
	}

	return qn;
}
