#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "util.h"

int* getclosedclasses(qnmodel* qn, int* numClosed)
{
	int r, count = 0;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] >= 0)
			count++;
	}
	*numClosed = count;
	if (count == 0) return NULL;
	int* idx = (int*)malloc(count * sizeof(int));
	int j = 0;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] >= 0)
			idx[j++] = r;
	}
	return idx;
}

int* getopenclasses(qnmodel* qn, int* numOpen)
{
	int r, count = 0;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] < 0)
			count++;
	}
	*numOpen = count;
	if (count == 0) return NULL;
	int* idx = (int*)malloc(count * sizeof(int));
	int j = 0;
	for (r = 0; r < qn->R; r++) {
		if (qn->N[r] < 0)
			idx[j++] = r;
	}
	return idx;
}
