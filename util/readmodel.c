#include <stdio.h>
#include <string.h>
#include "gmpla.h"
#include "util.h"
#include "assert.h"

qnmodel* readmodel(char* filename)
{
	qnmodel* qn=(qnmodel*)malloc(sizeof(qnmodel));
	int m;
	int r;

	/* initialize optional fields */
	qn->hasOpen = 0;
	qn->lambda = NULL;
	qn->isLD = 0;
	qn->Nt = 0;
	qn->mu = NULL;

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

	/* compute total closed population */
	qn->Nt = 0;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] > 0)
			qn->Nt += qn->N[r];
	}

	/* read optional keyword sections */
	char keyword[64];
	while (fscanf(f, "%63s", keyword) == 1) {
		if (strcmp(keyword, "LAMBDA") == 0) {
			qn->hasOpen = 1;
			qn->lambda = (mpq_t*)malloc(qn->R * sizeof(mpq_t));
			for (r = 0; r < qn->R; r++) {
				mpq_init(qn->lambda[r]);
				if (gmp_fscanf(f, "%Qd", &qn->lambda[r]) != 1) {
					printf("Error reading lambda[%d] from file\n", r);
					exit(1);
				}
				mpq_canonicalize(qn->lambda[r]);
			}
		} else if (strcmp(keyword, "MU") == 0) {
			qn->isLD = 1;
			if (qn->Nt <= 0) {
				printf("Error: MU section requires Nt > 0 (total closed population)\n");
				exit(1);
			}
			qn->mu = (mpq_t**)malloc(qn->M * sizeof(mpq_t*));
			for (m = 0; m < qn->M; m++) {
				qn->mu[m] = (mpq_t*)malloc(qn->Nt * sizeof(mpq_t));
				for (int k = 0; k < qn->Nt; k++) {
					mpq_init(qn->mu[m][k]);
					if (gmp_fscanf(f, "%Qd", &qn->mu[m][k]) != 1) {
						printf("Error reading mu[%d][%d] from file\n", m, k);
						exit(1);
					}
					mpq_canonicalize(qn->mu[m][k]);
				}
			}
		}
	}

	fclose(f);
	return qn;
}
