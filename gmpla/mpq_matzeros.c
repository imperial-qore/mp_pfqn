#include "gmpla.h"

mpq_t** mpq_matzeros(int rows, int cols)
{
    int i,j;
    mpq_t** matrix; 
    matrix = (mpq_t**) calloc(rows+1,sizeof(mpq_t)); 
        for(i = 0; i < rows+1; i++) 
        matrix[i] = (mpq_t*) calloc(cols+1,sizeof(mpq_t)); 
        
    for(i = 0; i < rows+1; i++) 
       for(j = 0; j < cols+1; j++)         
                mpq_init(matrix[i][j]);
    return matrix;
}
