/**********************************************************************
 *
 * File: gmdlj.c
 *
 * Create random ("gas") initial configuration for NVE Molecular
 * Dynamics of Lennard-Jonesium
 *
 * 04-May-2010 (MN)
 * 04-May-2012
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
  int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
  int i,m;
  double dr,dt,rho;
  double c;
  double *vx,*vy,*vz,*x,*y,*z;

  /* User input */

  for(;;) {
    printf("               n=");
    getval_i(&n);
    m=0;
    while(m*m*m<2*n) {   /* Check if magic number */
      m=m+2;
      if(m*m*m==2*n) {
	goto magic;
      }
    }
  }
  magic:

  printf("             rho=");
  getval_d(&rho);
  printf("              dt=");
  getval_d(&dt);
  printf("              dr=");
  getval_d(&dr);
  printf("         ntaskip=");
  getval_i(&ntaskip);
  printf(" ntprint/ntaskip=");
  getval_i(&ntprint);
  printf("         ntcskip=");
  getval_i(&ntcskip);
  printf("    ncor/ntcskip=");
  getval_i(&ncor);
  printf("  ntorig/ntcskip=");
  getval_i(&ntorig);
  printf("           ntjob=");
  getval_i(&ntjob);

  strcpy(fname,"mdlj_in.dat");
  printf("           fname=[mdlj_in.dat] ");
  getval_s(fname);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));
  vx=(double*)malloc(n*sizeof(double));
  vy=(double*)malloc(n*sizeof(double));
  vz=(double*)malloc(n*sizeof(double));

  /* Random positions & zero velocities */

  c=cbrt(n/rho);
  for(i=0;i<=n-1;i++) {
    x[i]=c*(drand48()-0.5);
    y[i]=c*(drand48()-0.5);
    z[i]=c*(drand48()-0.5);
    vx[i]=0.0;
    vy[i]=0.0;
    vz[i]=0.0;
  }

  /* Write startup file */

  nt=0;

  fpo=fopen(fname,"w");

  fwrite(&n,sizeof(int),1,fpo);
  fwrite(&ncor,sizeof(int),1,fpo);
  fwrite(&nt,sizeof(int),1,fpo);
  fwrite(&ntaskip,sizeof(int),1,fpo);
  fwrite(&ntcskip,sizeof(int),1,fpo);
  fwrite(&ntjob,sizeof(int),1,fpo);
  fwrite(&ntorig,sizeof(int),1,fpo);
  fwrite(&ntprint,sizeof(int),1,fpo);

  fwrite(&dr,sizeof(double),1,fpo);
  fwrite(&dt,sizeof(double),1,fpo);
  fwrite(&rho,sizeof(double),1,fpo);

  fwrite(x,sizeof(double),n,fpo);
  fwrite(y,sizeof(double),n,fpo);
  fwrite(z,sizeof(double),n,fpo);

  fwrite(vx,sizeof(double),n,fpo);
  fwrite(vy,sizeof(double),n,fpo);
  fwrite(vz,sizeof(double),n,fpo);

  fclose(fpo);

  /* Deallocate arrays */

  free(vx);
  free(vy);
  free(vz);
  free(x);
  free(y);
  free(z);

  return 0;

}

/**********************************************************************/
