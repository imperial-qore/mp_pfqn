#include "gmpla.h"

void mpq_vecdup(mpq_vec_t v1, mpq_vec_t v2, int n)
{
	int t;
	for (t=0;t<n;t++)
		mpq_set(v1[t],v2[t]);
}

