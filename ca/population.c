#include <stdio.h>
#include <stdlib.h>
#define bool int
#define true 1 
#define false 0
//#define DEBUG

int* initpop(int R)
{
	int* pop = (int*) calloc(R,sizeof(int));
	if (pop == NULL && R > 0) {
		fprintf(stderr, "Error: Failed to allocate memory for population array (%d elements)\n", R);
	}
	return pop;
}

int* resetpop(int* n,int R)
{
	free(n);
	return (int*) calloc(R,sizeof(int));
}

int nextpop(int* n, int* N, int R)
{
	int s=R;
	while (s>0 && (n[s-1]==N[s-1]))
		n[(s--)-1]=0;
	if (!s) return -1; // n equals N
	n[(s++)-1]++;
	return 0;
}

int nextpopinteractive(int* n, int* N, int R)
{
/*	int r;
	printf("current population: ");
	printpop(n);
	printf(" / ");
	printpop(N);
	for (r=1; r<=R;r++)
	{
		printf("N[%d]:=",r);
		scanf("%d",&v[r-1]);
	}
*/
	return 0;
}
int popindex(int* n,int R,int* planesizes) 
{
	int r,index=n[R-1];
	for (r=R-1;r>=1;r--)
		index += (n[r-1])*planesizes[r+1-1];
		
	return index;
}

int* getplanesizes(int *N, int R)
{
  int t;
  int* planesizes=(int *) calloc(R,sizeof(int));
  if (planesizes == NULL && R > 0) {
    fprintf(stderr, "Error: Failed to allocate memory for planesizes array (%d elements)\n", R);
    return NULL;
  }
  planesizes[R-1]=N[R-1]+1;
  for (t=R-1;t>=1;t--)
  	planesizes[t-1]=planesizes[(t+1)-1]*(N[t-1]+1);
  return planesizes;
}

void printpop(int* n,int R)
{
	int r;
	for(r=1;r<=R;r++)
		printf("%d ",n[r-1]);
}

bool issingleclass(int* n,int R)
{
	int nnz=0,r;
	for(r=1;r<=R;r++)
		if (n[r-1]>0)
			nnz++;
	if (nnz==1) return true;
	return false;
}

#ifdef DEBUG
int main()
{
	int R=3;
	int r;
	int *n=initpop(R);
	int N[]={1,1,2,'\0'};
	
	do	
	{
		printf("%d %d %d\n",n[0],n[1],n[2]);
	}
	while(!nextpop(n,N,R));

	return 0;
}
#endif
