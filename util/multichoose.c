#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"
#include "gmpla.h"

int ** multichoose(int n,int k)
{

    int **v=(int**)int_mat(nck(n+k-1,k),n,0); 
    if(v == NULL) {
        return NULL;
    }
    if(n==1)
    {
	v[0][0]=k;
	return v;
    }
    else if(k==0)
    {
    	return v;
    }
    else
    {

    	int i,last=0;
	for(i=0;i<=k;i++)
	{
		int t,j;
		long int nc=nck((n-1)+(k-i)-1,k-i);
		/* assign i choices to element 0 */
		for (t=0;t<nc;t++)
			v[last+t][0]=i;
		/* assign k-i choices to the rest */
    		int** w=multichoose(n-1,k-i);
		if (w == NULL) {
			free(v);
			return NULL;
		}
		for (t=0;t<nc;t++)
			for(j=1;j<=n-1;j++)
			{
				v[last+t][j]=w[t][j-1];
			}
		last+=nc;
		free(w);
	}
    	return v;
    }
}


