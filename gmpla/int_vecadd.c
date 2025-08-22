#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>

/* component-wise sum of two vectors  */
int* int_vecadd(int* v1, int* v2, int n)
{
        int* r=int_vec(n,0);
        for (;n>=1;n--)
                r[n-1]=v1[n-1]+v2[n-1];
        return r;
}

