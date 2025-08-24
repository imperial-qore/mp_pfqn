#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <gmpla.h>
#include "comom.h"

int hash(combsrep* Dn, int *comb, int i)
{
	if (i==1)
		return Dn->card*Dn->k+1+int_matmatchrow(Dn->combs,Dn->card,Dn->n, comb);
	else
		return (1+int_matmatchrow(Dn->combs,Dn->card,Dn->n,comb)-1)*Dn->k+i-1;
}


