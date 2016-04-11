/**********************************************************************
 *
 * File: gmclj.c
 *
 * Create random ("gas") initial configuration for NVT-Monte Carlo of
 * Lennard-Jonesium
 *
 * 12-May-2007 (MN)
 * 17-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.h"

#define NL '\n'
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
  double disp,dr,rho,t;
  double c;
  double *x,*y,*z;

  /* User input */

  printf("              n=");
  getval_i(&n);
  printf("            rho=");
  getval_d(&rho);
  printf("              t=");
  getval_d(&t);
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

  strcpy(fname,"mclj_in.dat");
  printf("          fname=[mclj_in.dat] ",fname);
  getval_s(fname);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Random positions */

  c=cbrt(n/rho);
  for(i=0;i<=n-1;i++) {
    x[i]=c*(drand48()-0.5);
    y[i]=c*(drand48()-0.5);
    z[i]=c*(drand48()-0.5);
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
  fwrite(&rho,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);

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
