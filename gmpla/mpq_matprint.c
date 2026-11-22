#include "gmpla.h"

void mpq_matprint(mpq_mat_t mat, int n, int m)
{
	int i,j;
	for (i=1;i<=n;i++)
	{
		for (j=1;j<=m;j++)
		{
			printf("%3g ",(double)mpq_get_d(mat[i-1][j-1]));
		}
		printf("\n");
	}
	printf("\n");
}


