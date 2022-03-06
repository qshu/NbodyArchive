/*
 * mysleep fortran interface to sleep(3)
 */

#include <sys/types.h>
#include <sys/times.h>
#ifdef ibm
#define mysleep_ mysleep
#endif

void mysleep_(time)
    int * time;
{
    sleep(*time);
}
