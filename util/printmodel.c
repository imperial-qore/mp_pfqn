#include <stdio.h>
#include "util.h"

void printmodel(qnmodel* qn)
{
	int m,r,mtot=0;
	for(m=1;m<=qn->M;m++)
		mtot+=qn->mi[m-1];
	printf("The queueing network model has %d queues (%d are replicas) and %d classes\n",mtot,mtot-qn->M,qn->R);
	printf("\t  N[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++)
		printf("%12.5f",(double)qn->N[r-1]);
	printf("\n");
	printf("\t  Z[1:%d]:",qn->R);
	for (r=1;r<=qn->R;r++)
		printf("%12.5f",(double)qn->Z[r-1]);
	printf("\n");
	for(m=1;m<=qn->M;m++)
	{
	 printf("mi=%d    L[%d,1:%d]:",qn->mi[m-1],m,qn->R);
	 for (r=1;r<=qn->R;r++)
	 {
		printf("%12.5f",(double)qn->L[m-1][r-1]);
	 }
	 printf("\n");
	}
	printf("\n");
}

/*
int main(int argc, char** argv)
{
	printf("Opening %s\n",argv[1]);
	qnmodel* qn=(qnmodel*) readmodel(argv[1]);
	printmodel(qn);
	return 0;
}
*/
