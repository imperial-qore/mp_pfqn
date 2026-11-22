#include "gmp.h"

int nnz(int* v, int n)
{
	int i;
	int nz=0; /* number of non-zeros */
	for (i=0;i<n;i++)
	{
		if (v[i]!=0) nz++;
	}
	return nz;
}


