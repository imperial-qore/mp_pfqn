#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "util.h" 

int randi(int min,int max)
{
	return min + (rand() % (max-min));
}


