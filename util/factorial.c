#include "util.h"

long int factorial(long int N)
{
	if (N<0) return -1;
	if (N==0) return 1;
	if (N==1) return 1;
	return (long int) N*factorial(N-1);
}
