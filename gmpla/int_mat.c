#include "gmpla.h" 

int** int_mat(int rows, int cols, int initval)
{
    int i; 
    int j; 
    int** matrix; 
    matrix = (int**) calloc(rows,sizeof(int*)); 
    if (matrix == NULL) {
        return NULL;
    }
    for(i = 0; i < rows; i++) {
        matrix[i] = (int*) calloc(cols,sizeof(int)); 
        if (matrix[i] == NULL) {
            // Clean up previously allocated memory
            for (j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
        
    for(i = 0; i < rows; i++) 
        for(j = 0; j < cols; j++)        { 
                matrix[i][j]=initval;
	}
    return matrix;
}
