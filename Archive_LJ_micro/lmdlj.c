/**********************************************************************
 *
 * File: lmdlj.c
 *
 * Create initial configuration (fcc lattice) for NVE Molecular Dynamics
 * of Lennard-Jonesium
 * Random velocities
 *
 * 08-May-2010 (MN)
 * 04-May-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.c"

#define BSIZE 80

/**********************************************************************/

int main() {

/**********************************************************************/

  FILE *fpo;
  char fname[BSIZE];
  int n,ncor,nt,ntaskip,ntcskip,ntjob,ntprint,ntorig;
  int i,ix,iy,iz,m;
  double dr,dt,rho;
  double c,fact,sv2,svx,svy,svz,t;
  double *vx,*vy,*vz,*x,*y,*z;

  /* User input */

  for(;;) {
    printf("               n=");
    getval_i(&n);
    m=0;
    while(m*m*m<2*n) {
      m=m+2;
      if(m*m*m==2*n) {   /* Check if magic number */
	goto magic;
      }
    }
  }
  magic:

  printf("             rho=");
  getval_d(&rho);
  printf("               t=");
  getval_d(&t);
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

  vx=(double*)malloc(n*sizeof(double));
  vy=(double*)malloc(n*sizeof(double));
  vz=(double*)malloc(n*sizeof(double));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Fcc lattice & random velocities */

  c=cbrt(n/rho);
  i=0;
  for(ix=0;ix<=m-1;ix++) {
    for(iy=0;iy<=m-1;iy++) {
      for(iz=0;iz<=m-1;iz++) {
	if((ix+iy+iz)%2==0) {
	  x[i]=c*((ix+0.5)/m-0.5);
	  y[i]=c*((iy+0.5)/m-0.5);
	  z[i]=c*((iz+0.5)/m-0.5);
	  vx[i]=drand48()-0.5;
	  vy[i]=drand48()-0.5;
	  vz[i]=drand48()-0.5;
	  i=i+1;
	}
      }
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
