/**********************************************************************
 *
 * File: ggemclj.c
 *
 * Create random ("gas") initial configuration for Gibbs Ensemble (NVT)
 * Monte Carlo of Lennard-Jonesium
 *
 * 01-May-2010 (MN)
 * 19-Apr-2012
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
  int n,n1,nt,ntjob,ntprint,ntskip;
  int i,n2;
  double disp1,disp2,dv,ptr,pvm,t,v,v1;
  double c1,c2,v2;
  double *x,*y,*z;

  /* User input */

  printf("             n1=");
  getval_i(&n1);
  printf("             n2=");
  getval_i(&n2);
  printf("             v1=");
  getval_d(&v1);
  printf("             v2=");
  getval_d(&v2);
  printf("              t=");
  getval_d(&t);
  printf("          disp1=");
  getval_d(&disp1);
  printf("          disp2=");
  getval_d(&disp2);
  printf("             dv=");
  getval_d(&dv);

  for(;;) {
    printf("    prob(vmove)=");
    getval_d(&pvm);
    printf("    prob(trans)=");
    getval_d(&ptr);
    if(ptr+pvm<1.0) {
      break;
    }
  }

  printf("         ntskip=");
  getval_i(&ntskip);
  printf(" ntprint/ntskip=");
  getval_i(&ntprint);
  printf("   ntjob/ntskip=");
  getval_i(&ntjob);

  strcpy(fname,"gemclj_in.dat");
  printf("          fname=[gemclj_in.dat] ");
  getval_s(fname);

  /* Parameters */

  n=n1+n2;
  c1=cbrt(v1);
  c2=cbrt(v2);
  v=v1+v2;

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Random positions: Box 1 */

  for(i=0;i<=n1-1;i++) {
    x[i]=(drand48()-0.5)*c1;
    y[i]=(drand48()-0.5)*c1;
    z[i]=(drand48()-0.5)*c1;
  }

  /* Box 2 */

  for(i=n1;i<=n-1;i++) {
    x[i]=(drand48()-0.5)*c2;
    y[i]=(drand48()-0.5)*c2;
    z[i]=(drand48()-0.5)*c2;
  }

  /* RNG seed */

  seed=seed48(dummy);

  /* Write startup file */

  nt=0;

  fpo=fopen(fname,"w");

  fwrite(&n,sizeof(int),1,fpo);
  fwrite(&n1,sizeof(int),1,fpo);
  fwrite(&nt,sizeof(int),1,fpo);
  fwrite(&ntjob,sizeof(int),1,fpo);
  fwrite(&ntprint,sizeof(int),1,fpo);
  fwrite(&ntskip,sizeof(int),1,fpo);
  fwrite(seed,sizeof(unsigned short),3,fpo);

  fwrite(&disp1,sizeof(double),1,fpo);
  fwrite(&disp2,sizeof(double),1,fpo);
  fwrite(&dv,sizeof(double),1,fpo);
  fwrite(&ptr,sizeof(double),1,fpo);
  fwrite(&pvm,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);
  fwrite(&v,sizeof(double),1,fpo);
  fwrite(&v1,sizeof(double),1,fpo);

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
