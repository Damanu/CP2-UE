/**********************************************************************
 *
 * File: gnptmclj.c
 *
 * Create random ("gas") initial configuration for NpT-Monte Carlo of
 * Lennard-Jonesium
 *
 * 03-Apr-2010 (MN)
 * 17-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.c"

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
  int i,m;
  double disp,dlnv,dr,p,pvm,t,v;
  double c,rho;
  double *x,*y,*z;

  /* User input */

  for(;;) {
    printf("              n=");
    getval_i(&n);
    m=0;
    while(m*m*m<2*n) {
      m=m+2;
      if(m*m*m==2*n) { /* Check if magic number */
        goto magic;
      }
    }
  }
  magic:

  printf("              p=");
  getval_d(&p);
  printf("              t=");
  getval_d(&t);
  printf("            rho=");
  getval_d(&rho);
  printf("           disp=");
  getval_d(&disp);
  printf("           dv/v=");
  getval_d(&dlnv);
  printf("    prob(vmove)=");
  getval_d(&pvm);
  printf("             dr=");
  getval_d(&dr);
  printf("         ntskip=");
  getval_i(&ntskip);
  printf(" ntprint/ntskip=");
  getval_i(&ntprint);
  printf("   ntjob/ntskip=");
  getval_i(&ntjob);

  strcpy(fname,"nptmclj_in.dat");
  printf("          fname=[nptmclj_i.dat] ");
  getval_s(fname);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Random positions */

  v=n/rho;
  c=cbrt(v);
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
  fwrite(&seed,sizeof(unsigned short),3,fpo);

  fwrite(&disp,sizeof(double),1,fpo);
  fwrite(&dlnv,sizeof(double),1,fpo);
  fwrite(&dr,sizeof(double),1,fpo);
  fwrite(&p,sizeof(double),1,fpo);
  fwrite(&pvm,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);
  fwrite(&v,sizeof(double),1,fpo);

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
