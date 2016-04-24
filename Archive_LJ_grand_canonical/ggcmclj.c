/**********************************************************************
 *
 * File: ggcmclj.c
 *
 * Create random ("gas") initial configuration for Grand Canonical
 * (TVmu) Monte Carlo of Lennard-Jonesium
 *
 * 08-Apr-2010 (MN)
 * 18-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.h"

#define BSIZE 80


/**********************************************************************/

int main() {

/**********************************************************************/

  FILE *fpo;
  char fname[BSIZE];
  unsigned short dummy[3];
  unsigned short *seed;
  int n,nt,ntjob,ntprint,ntskip;
  int i;
  double pcre,disp,dr,t,v,zz;
  double c;
  double *x,*y,*z;

  /* User input */

  printf("              t=");
  getval_d(&t);
  printf("              v=");
  getval_d(&v);
  printf("              z=");
  getval_d(&zz);
  printf("              n=");
  getval_i(&n);

  for(;;) {
    printf("     p(cre/del)=");
    getval_d(&pcre);
    if(pcre<0.5) {
      break;
    }
  }

  printf("           disp=");
  getval_d(&disp);
  printf("             dr=");
  getval_d(&dr);
  printf("         ntskip=");
  getval_i(&ntskip);
  printf(" ntprint/ntskip=");
  getval_i(&ntprint);
  printf("   ntjob/ntskip=");
  getval_i(&ntjob);

  strcpy(fname,"gcmclj_in.dat");
  printf("          fname=[gcmclj_in.dat] ");
  getval_s(fname);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Box length */

  c=cbrt(v);

  /* Random particle positions */

  for(i=0;i<=n-1;i++) {
    x[i]=(drand48()-0.5)*c;
    y[i]=(drand48()-0.5)*c;
    z[i]=(drand48()-0.5)*c;
  }

  /* RNG seed */

  seed=seed48(dummy);

  /* Write startup file */

  nt=0;

  fpo=fopen(fname,"w");

  fwrite(&n,sizeof(int),1,fpo);
  fwrite(&nt,sizeof(int),1,fpo);
  fwrite(&ntjob,sizeof(int),1,fpo);
  fwrite(&ntprint,sizeof(int),1,fpo);
  fwrite(&ntskip,sizeof(int),1,fpo);
  fwrite(seed,sizeof(unsigned short),3,fpo);

  fwrite(&disp,sizeof(double),1,fpo);
  fwrite(&dr,sizeof(double),1,fpo);
  fwrite(&pcre,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);
  fwrite(&v,sizeof(double),1,fpo);
  fwrite(&zz,sizeof(double),1,fpo);

  fwrite(x,sizeof(double),n,fpo);
  fwrite(y,sizeof(double),n,fpo);
  fwrite(z,sizeof(double),n,fpo);

  fclose(fpo);

  /* Deallocate arrays */

  free(x);
  free(y);
  free(z);

  return 0;

}

/**********************************************************************/
