#include "gmp.h"

int nnzposcmp(int*i1, int*i2, int n)
{
	int nnz1=nnz(i1,n);
	int nnz2=nnz(i2,n);

	if (nnz1>nnz2)
	{
		return 1;
	}
	else if(nnz1<nnz2)
	{
		return 0;
	}
	else /* nnz1==nnz2 */
	{
		int j;
		for (j=0;j<n;j++)	
		{
			if (i1[j]==0 && i2[j]>0)	
			{
				return 1;		
			}
			else if (i1[j]>0 && i2[j]==0)
			{
				return 0;
			}
		}
		return 0;
	}
}


