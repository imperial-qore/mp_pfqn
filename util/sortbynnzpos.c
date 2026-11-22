#include "gmp.h"
#include "util.h" 

int** sortbynnzpos(int** I, int m,int n)
{
	int i,j;
	int *t;
	for (i=0;i<m-1;i++)
	{
		for (j=i+1;j<m;j++)
		{
			if (nnzposcmp(I[i],I[j],n)==1)
			{
				t=I[i];
				I[i]=I[j];
				I[j]=t;
			}
		}
	}
	return I;
}


