// Main file with functions

#include "membranes/membrane.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef INT_MAX
#include <limits.h>
#endif

void err_exit(int code, const char *fmt, ...)
{
    va_list argp;
    FILE *fp;
    char msg[1024];

    va_start(argp, fmt);
    vsnprintf(msg, 1024, fmt, argp);
    va_end(argp);

    fp = stdout;
    fprintf(fp, "\nERROR EXIT with code %d ==> %s\n\n", code, msg);
    fflush(fp);

    exit(code);
}

void err_locate(int l, const char *f, const char *d, const char *t, const char *u)
{}

#define MAX_WARNING_CODE 100
#define MAX_WARNINGS     30

int warn_cnts[MAX_WARNING_CODE+1] = { [0] = 0 };

void Warning(int code, const char *fmt, ...)
{
  va_list argp;
  char    buf[256];

  if( code > MAX_WARNING_CODE || code < 0 ) {
    Warning( 1, "warning code %d is out of range; changed to 1", code );
    code = 1;
  }

  if( warn_cnts[code] < INT_MAX )
    ++warn_cnts[code];

  if( warn_cnts[code] > MAX_WARNINGS )
    return;

  if( warn_cnts[code] == MAX_WARNINGS )
    snprintf( buf, sizeof(buf), "\nLast warning from process %d: %s (code %d)\n", 
                                0, fmt, code );
  else
    snprintf( buf, sizeof(buf), "\nWarning on process %d: %s (code %d)\n\n", 
                                0, fmt, code );

  va_start  ( argp, fmt );
  vfprintf  ( stdout, buf, argp );
  va_end    ( argp );
} 

char *node_to_string(long n)
{
   char *s = malloc(1024), *r;
   snprintf(s, 1024, "node %ld vertex", n);
   r = strdup(s);
   free(s);
   return r;
}

void ion_infinite(Membrane_info *ifo,
                  Membrane_cell_info *ci,
                  vm_t* Vmvec,
                  yyy_t* yyvec,
                  const int nnodes)
{
    int nsvar = ifo->Nsvar;
    for (int n=0; n<nnodes; ++n)
    {
        (*ifo->infinite)(ci, &Vmvec[n], &yyvec[n * nsvar], 1);
    }
}

void ion_step(Membrane_info *ifo,
              Membrane_cell_info *ci,
              vm_t   *Vmvec,
              yyy_t  *yyvec,
              vm_t   *Isd,
              float  *dtime,
              double *Imi,
              const double dt,
              const double simtime,
              const int nnodes)
{
    int nsvar = ifo->Nsvar;
    
    #pragma omp parallel for schedule(static)
    for (int n=0; n<nnodes; ++n)
    {
        Imi[n] = (*ifo->step)(Vmvec[n], &yyvec[n * nsvar], ci, Isd[n],
                              dt, &dtime[n], simtime, n);
    }
}
