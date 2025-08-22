#include "gmpla.h"

void mpq_vecprint(mpq_vec_t vec, int n)
{
	int i;
	for (i=1;i<=n;i++)
	{
		gmp_printf("%Q12d ",vec[i-1]);
		printf("\n");
	}
	printf("\n");
}

