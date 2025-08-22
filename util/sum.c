#include "util.h" 

int sum(int* v, int n)
{
	int i;
	int S=0;
	for(i=0;i<n;i++)
		S+=v[i];

	return S;
}


