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

	/* sanity checks on the fixed part of the model */
	if (qn->R < 1 || qn->M < 1) {
		printf("Error: model needs R >= 1 and M >= 1 (got R=%d, M=%d)\n", qn->R, qn->M);
		exit(1);
	}
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] < -1) {
			printf("Error: N[%d]=%d is invalid (use -1 for an open class)\n", r, qn->N[r]);
			exit(1);
		}
		if (mpz_sgn(qn->Z[r]) < 0) {
			printf("Error: think time Z[%d] is negative\n", r);
			exit(1);
		}
	}
	for (m = 0; m < qn->M; m++) {
		if (qn->mi[m] < 1) {
			printf("Error: multiplicity mi[%d]=%d must be >= 1\n", m, qn->mi[m]);
			exit(1);
		}
		for (r = 0; r < qn->R; r++) {
			if (mpz_sgn(qn->L[m][r]) < 0) {
				printf("Error: service demand L[%d][%d] is negative\n", m, r);
				exit(1);
			}
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
				if (mpz_sgn(mpq_denref(qn->lambda[r])) == 0) {
					printf("Error: lambda[%d] has a zero denominator\n", r);
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
					if (mpz_sgn(mpq_denref(qn->mu[m][k])) == 0) {
						printf("Error: mu[%d][%d] has a zero denominator\n", m, k);
						exit(1);
					}
					mpq_canonicalize(qn->mu[m][k]);
				}
			}
		} else {
			/* Anything else here is a malformed file: an extra station
			 * line, an extra column on a station line, or a misspelled
			 * section keyword.  Silently skipping it parses a different
			 * model from the one the user wrote. */
			printf("Error: unexpected token '%s' after the station block (check M, R and the LAMBDA/MU keywords)\n", keyword);
			exit(1);
		}
	}

	fclose(f);
	return qn;
}
