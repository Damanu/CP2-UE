/**********************************************************************
 *
 * File: zmdlj.c
 *
 * NVE Molecular Dynamics of Lennard-Jonesium
 * Re-initialize checkpoint file
 *
 * 07-May-2010 (MN)
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

  FILE *fpi,*fpo;
  char fname[BSIZE];
  int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
  int i;
  double dr,dt,rho;
  double fact,rhon,sk,skn,svx,svy,svz;
  double *vx,*vy,*vz,*x,*y,*z;

  /* Read checkpoint file */

  strcpy(fname,"mdlj_out.dat");
  printf("          infile=[mdlj_out.dat] ");
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

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));
  vx=(double*)malloc(n*sizeof(double));
  vy=(double*)malloc(n*sizeof(double));
  vz=(double*)malloc(n*sizeof(double));

  /* Positions & velocities */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(vx,sizeof(double),n,fpi);
  fread(vy,sizeof(double),n,fpi);
  fread(vz,sizeof(double),n,fpi);

  fclose(fpi);

  /* Save old density & kinetic energy */

  rhon=rho;
  sk=0.0;
  for(i=0;i<=n-1;i++) {
    sk=sk+(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }
  sk=0.5*sk/n;
  skn=sk;

  /* User input */

  printf("               n=%13d\n",n);
  printf("             rho=[%12.5lf] ",rho);
  getval_d(&rhon);
  printf("             k/n=[%12.5lf] ",sk);
  getval_d(&skn);
  printf("              dt=[%12.5lf] ",dt);
  getval_d(&dt);
  printf("              dr=[%12.5lf] ",dr);
  getval_d(&dr);
  printf("         ntaskip=[%12d] ",ntaskip);
  getval_i(&ntaskip);
  printf(" ntprint/ntaskip=[%12d] ",ntprint);
  getval_i(&ntprint);
  printf("         ntcskip=[%12d] ",ntcskip);
  getval_i(&ntcskip);
  printf("    ncor/ntcskip=[%12d] ",ncor);
  getval_i(&ncor);
  printf("  ntorig/ntcskip=[%12d] ",ntorig);
  getval_i(&ntorig);
  printf("           ntjob=[%12d] ",ntjob);
  getval_i(&ntjob);

  /* Rescale positions */

  if(rhon!=rho) {
    fact=cbrt(rho/rhon);
    rho=rhon;
    for(i=0;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }
  }

  /* Zero net momentum & rescale velocities */

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
  for(i=0;i<=n-1;i++) {
    vx[i]=vx[i]-svx;
    vy[i]=vy[i]-svy;
    vz[i]=vz[i]-svz;
  }

  if(skn!=sk) {
    fact=sqrt(skn/sk);
    for(i=0;i<=n-1;i++) {
      vx[i]=fact*vx[i];
      vy[i]=fact*vy[i];
      vz[i]=fact*vz[i];
    }
  }

  /* Write new startup file */

  strcpy(fname,"mdlj_in.dat");
  printf("         outfile=[ mdlj_in.dat] ");
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

