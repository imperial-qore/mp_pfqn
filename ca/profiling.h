#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
#define MEMUSAGE (getrusage(RUSAGE_CHILDREN,&ruse), ruse.ru_maxrss)
//#define CPUTIME 1 
int AORSCTR; // addition or subtraction counter
int MULCTR; // multiplcation counter
int DIVCTR; // division counter
double AORSTIME; // addition or subtraction time counter
double MULTIME; // multiplcation time counter
double DIVTIME; // division time counter
double t0,t1;
struct rusage ruse;
extern int getrusage();
