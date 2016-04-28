/**********************************************************************
 *
 * File: zisomdlj.c
 *
 * Gaussian isokinetic (NVT) Molecular Dynamics of Lennard-Jonesium
 * Re-initialize checkpoint file
 *
 * 15-Nov-2007 (MN)
 * 05-May-2012
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

  FILE *fpi,*fpo;
  char fname[BSIZE];
  int i;
  int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
  double fact,rhon,sv2,svx,svy,svz;
  double dr,dt,rho,t;
  double *vx,*vy,*vz,*x,*y,*z;

  /* Read checkpoint file */

  strcpy(fname,"isomdlj_out.dat");
  printf("          infile=[isomdlj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  /* Simulation parameters */

  fread(&n,sizeof(int),1,fpi);
  fread(&ncor,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntaskip,sizeof(int),1,fpi);
  fread(&ntcskip,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntorig,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);

  fread(&dr,sizeof(double),1,fpi);
  fread(&dt,sizeof(double),1,fpi);
  fread(&rho,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);

  /* Allocate arrays */

  vx=(double*)malloc(n*sizeof(double));
  vy=(double*)malloc(n*sizeof(double));
  vz=(double*)malloc(n*sizeof(double));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions & velocities */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(vx,sizeof(double),n,fpi);
  fread(vy,sizeof(double),n,fpi);
  fread(vz,sizeof(double),n,fpi);

  fclose(fpi);

  rhon=rho;

  /* User input */

  printf("               n=%16d\n",n);
  printf("             rho=[%15.5lf] ",rho);
  getval_d(&rhon);
  printf("               t=[%15.5lf] ",t);
  getval_d(&t);
  printf("              dt=[%15.5lf] ",dt);
  getval_d(&dt);
  printf("              dr=[%15.5lf] ",dr);
  getval_d(&dr);
  printf("         ntaskip=[%15d] ",ntaskip);
  getval_i(&ntaskip);
  printf(" ntprint/ntaskip=[%15d] ",ntprint);
  getval_i(&ntprint);
  printf("         ntcskip=[%15d] ",ntcskip);
  getval_i(&ntcskip);
  printf("    ncor/ntcskip=[%15d] ",ncor);
  getval_i(&ncor);
  printf("  ntorig/ntcskip=[%15d] ",ntorig);
  getval_i(&ntorig);
  printf("           ntjob=[%15d] ",ntjob);
  getval_i(&ntjob);

  /* Rescale positions */

  if(rhon/=rho) {
    fact=cbrt(rho/rhon);
    for(i=0;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }
  }

  /* Zero net momentum & set initial temperature */

  svx=0.0;
  svy=0.0;
  svz=0.0;
  for(i=0;i<=n-1;i++) {
    svx=svx+vx[i];
    svy=svy+vy[i];
    svz=svz+vz[i];
  }
  svx=svx/n;
  svy=svy/n;
  svz=svz/n;

  sv2=0.0;
  for(i=0;i<=n-1;i++) {
    vx[i]=vx[i]-svx;
    vy[i]=vy[i]-svy;
    vz[i]=vz[i]-svz;
    sv2=sv2+(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }

  fact=sqrt(3*n*t/sv2);
  for(i=0;i<=n-1;i++) {
    vx[i]=fact*vx[i];
    vy[i]=fact*vy[i];
    vz[i]=fact*vz[i];
  }

  /* Write new startup file */

  strcpy(fname,"isomdlj_in.dat");
  printf("         outfile=[ isomdlj_in.dat] ");
  getval_s(fname);

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
  fwrite(&t,sizeof(double),1,fpo);

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
