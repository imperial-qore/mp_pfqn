#include "gmpla.h"

void int_vecprint(int* vec, int n)
{
	int i;
	for (i=1;i<=n;i++)
		printf("%3d ",vec[i-1]);
	printf("\n");
}

