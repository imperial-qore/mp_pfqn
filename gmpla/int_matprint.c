#include "gmpla.h"

void int_matprint(int** mat, int n, int m)
{
	int i,j;
	for (i=0;i<n;i++)
	{
		for (j=0;j<m;j++)
			printf("%3d ",mat[i][j]);
		printf("\n");
	}
}


