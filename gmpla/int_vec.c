#include <stdlib.h>

int* int_vec(int n, int value)
{
    int* vec=calloc(n+1,sizeof(int));
    int i;
    for (i=0;i<n;i++)
     vec[i]=value;
    return vec;
}
