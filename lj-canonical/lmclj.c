/**********************************************************************
 *
 * File: lmclj.c
 *
 * Create initial configuration (fcc lattice) for NVT-Monte Carlo of
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
  unsigned short dummy[3],*seed;
  int n,nt,ntjob,ntprint,ntskip;
  int i,ix,iy,iz,m;
  double disp,dr,rho,t;
  double c;
  double *x,*y,*z;

  /* User input */

  for(;;) {
    printf("              n=");
    scanf("%d",&n);
    m=0;
    while(m*m*m<2*n) {
      m=m+2;
      if(m*m*m==2*n) { /* Check if magic number */
        goto magic;
      }
    }
  }
  magic:

  printf("            rho=");
  scanf("%lf",&rho);
  printf("              t=");
  scanf("%lf",&t);
  printf("           disp=");
  scanf("%lf",&disp);
  printf("             dr=");
  scanf("%lf",&dr);
  printf("         ntskip=");
  scanf("%d",&ntskip);
  printf(" ntprint/ntskip=");
  scanf("%d",&ntprint);
  printf("   ntjob/ntskip=");
  scanf("%d",&ntjob);

  strcpy(fname,"mclj_in.dat");

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Fcc lattice */

  c=cbrt(n/rho);
  i=0;
  for(ix=0;ix<=m-1;ix++) {
    for(iy=0;iy<=m-1;iy++) {
      for(iz=0;iz<=m-1;iz++) {
	if((ix+iy+iz)%2==0) {
 	  x[i]=c*((ix-0.5)/m-0.5);
 	  y[i]=c*((iy-0.5)/m-0.5);
 	  z[i]=c*((iz-0.5)/m-0.5);
	  i=i+1;
	}
      }
    }
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
