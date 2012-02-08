/* Concurrent read version */

/*  $Id: htime.c 852 2011-10-27 15:14:53Z wrp $ */
/* $Revision: 852 $  */

/* modify to use getrusage rather than time() or times() */

#include <stdio.h>
#include <time.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/time.h>
#ifdef TIMES
#include <sys/times.h>
#else
#undef TIMES
#endif
#endif

#ifndef HZ
#define HZ 100
#endif

long s_time ()			/* returns time in milliseconds */
{
#ifndef TIMES
  time_t time(), tt;		/* time() - returns time in seconds */
  return time(&tt)*1000;
#else
  struct tms tt;
  times(&tt);		/* times() returns tv.tms_utime in CLK_TCK's */
#ifdef CLK_TCK
  return tt.tms_utime*1000/CLK_TCK;
#else
  return tt.tms_utime*1000/HZ;
#endif
#endif
}

void ptime (FILE *fp, long time)		/* prints the time */
{
  fprintf (fp, "%6.3f",(double)(time)/1000.0);
}
