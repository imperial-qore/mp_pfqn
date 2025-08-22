#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

/* int_matmatchrow - find row in integer matrix 
   INPUT:
   mat   - integer matrix
   nrows - number of rows in mat
   ncols - number of cols in mat
   row   - target row vector
   OUTPUT: 
   -1    - no row found
   int   - row index of first match  
*/
long int int_matmatchrow(int** mat, long int nrows, long int ncols, int* row)
{
	long int i,j,identical;
	for(i=0;i<nrows;i++)
	{
		identical=1;
		for(j=0;j<ncols;j++)
			if (row[j]!=mat[i][j]) 
			{
				identical=0;
				break;
			}
		if (identical) return i;
	}
	return -1;
}
