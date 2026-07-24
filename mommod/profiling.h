#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
#define MEMUSAGE (getrusage(RUSAGE_CHILDREN,&ruse), ruse.ru_maxrss)
//#define CPUTIME 1 
extern int AORSCTR; // addition or subtraction counter
extern int MULCTR; // multiplcation counter
extern int DIVCTR; // division counter
extern double AORSTIME; // addition or subtraction time counter
extern double MULTIME; // multiplcation time counter
extern double DIVTIME; // division time counter
extern double t0,t1;
extern struct rusage ruse;
