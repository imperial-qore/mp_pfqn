/* Globals expected by profiling.h and the shared mommod sources.  Kept
 * apart from main.c so that the engine bundle linked into bin/mom carries
 * its own (local) copies. */
#include <stdbool.h>
#include <gmp.h>
#include <zpla.h>
#include "mommod.h"

/* globals expected by profiling.h and the shared solver sources */
int AORSCTR;
int MULCTR;
int DIVCTR;
double AORSTIME;
double MULTIME;
double DIVTIME;
double t0,t1;
struct rusage ruse;

qnmodel* qnm;
bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;
int VERBOSE_MOD = 0;
