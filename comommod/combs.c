#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"

/* The basis combinations for a class depend only on (r,M), never on the
 * prime, so they are built once in main and shared read-only by every
 * modular image.  comom rebuilds them per class and never frees them,
 * which is harmless for a single solve but would leak once per prime
 * here. */

/* Mixed-radix encoding of a combination, each coordinate in [0,M].
 * Returns 0 when the encoding cannot be represented, in which case the
 * lookup is not built and hash() falls back to the linear scan. */
static int code_fits(int n, int M)
{
	double e = 1.0;
	int j;
	for (j=0;j<n;j++)
	{
		e *= (double)(M+1);
		if (e > 9.0e18) return 0;
	}
	return 1;
}

static unsigned long long encode(int* comb, int n, int M)
{
	unsigned long long c = 0, w = 1;
	int j;
	for (j=0;j<n;j++)
	{
		c += (unsigned long long) comb[j] * w;
		w *= (unsigned long long)(M+1);
	}
	return c;
}

static int cmp_key(const void* a, const void* b)
{
	unsigned long long x = ((const unsigned long long*)a)[0];
	unsigned long long y = ((const unsigned long long*)b)[0];
	return (x<y) ? -1 : (x>y);
}

/* Build the class-r basis exactly as comom/main.c does: multichoose(r,M),
 * the last coordinate forced to zero, then sorted by nnz position. */
combsrep* makeDn(int r, int M)
{
	combsrep* Dn = (combsrep*) calloc(1,sizeof(combsrep));
	int i;

	if (Dn == NULL) return NULL;
	Dn->n = r;
	Dn->k = M;
	Dn->card = (int) nck(M+r-1,M);
	Dn->combs = (int**) multichoose(Dn->n,Dn->k);
	if (Dn->combs == NULL) { free(Dn); return NULL; }
	for (i=1;i<=Dn->card;i++)
		Dn->combs[i-1][r-1] = 0;
	Dn->combs = (int**) sortbynnzpos(Dn->combs,Dn->card,Dn->n);

	if (code_fits(Dn->n, M))
	{
		/* pairs of (code, row) sorted by code, then split */
		unsigned long long* pair = (unsigned long long*)
			malloc((size_t)Dn->card * 2 * sizeof(unsigned long long));
		if (pair != NULL)
		{
			for (i=0;i<Dn->card;i++)
			{
				pair[2*i]   = encode(Dn->combs[i], Dn->n, M);
				pair[2*i+1] = (unsigned long long) i;
			}
			qsort(pair, (size_t)Dn->card, 2*sizeof(unsigned long long), cmp_key);
			Dn->keys  = (unsigned long long*) malloc((size_t)Dn->card*sizeof(unsigned long long));
			Dn->order = (int*) malloc((size_t)Dn->card*sizeof(int));
			if (Dn->keys && Dn->order)
				for (i=0;i<Dn->card;i++)
				{
					Dn->keys[i]  = pair[2*i];
					Dn->order[i] = (int) pair[2*i+1];
				}
			else
			{
				free(Dn->keys); free(Dn->order);
				Dn->keys = NULL; Dn->order = NULL;
			}
			free(pair);
		}
	}
	return Dn;
}

void freeDn(combsrep* Dn)
{
	int i;
	if (!Dn) return;
	if (Dn->combs)
	{
		for (i=0;i<Dn->card;i++) free(Dn->combs[i]);
		free(Dn->combs);
	}
	free(Dn->keys);
	free(Dn->order);
	free(Dn);
}

/* hash - position of a basis entry, identical to comom/hash.c.
 *
 *   i == 0  : the G value for this population, stored after all G^k
 *   i in 1..M : the G^k value, stored at population_position*M + (i-1)
 *
 * The only difference from comom is that the population lookup uses the
 * precomputed ordering when it is available; the returned index is the
 * same either way. */
int hash(combsrep* Dn, int *comb, int i)
{
	long int population_position;

	if (Dn->keys != NULL)
	{
		unsigned long long c = encode(comb, Dn->n, Dn->k);
		int lo = 0, hi = Dn->card - 1;
		population_position = -1;
		while (lo <= hi)
		{
			int mid = lo + (hi-lo)/2;
			if (Dn->keys[mid] == c) { population_position = Dn->order[mid]; break; }
			if (Dn->keys[mid] < c) lo = mid+1; else hi = mid-1;
		}
	}
	else
		population_position = int_matmatchrow(Dn->combs, Dn->card, Dn->n, comb);

	if (population_position < 0)
	{
		fprintf(stderr, "Error: hash() - combination not found in basis\n");
		return -1;
	}

	if (i == 0)
		return Dn->card * Dn->k + (int) population_position + 1;
	return (int) population_position * Dn->k + (i - 1) + 1;
}
