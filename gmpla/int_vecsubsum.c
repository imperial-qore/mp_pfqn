#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>

/* sum of subvector elements */
int int_vecsubsum(int* v, int from, int to)
{
        int r=0;

        for (;from<=to;from++)
                r = r + v[from];

        return r;
}

