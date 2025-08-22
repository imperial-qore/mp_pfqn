#include "gmpla.h"

mpq_vec_t mpq_vec(int rows,int valuenum, int valueden)
{
    int i;
    mpq_t* vec; 
    vec = (mpq_t*) calloc(rows+1,sizeof(mpq_t)); 
        
    for(i=0;i<rows+1;i++) 
    {
	mpq_init(vec[i]);
	mpq_set_si(vec[i],valuenum,valueden);
    }
    return (mpq_vec_t) vec;
}
