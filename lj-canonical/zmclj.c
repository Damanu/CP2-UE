/**********************************************************************
 *
 * File: zmclj.c
 *
 * NVT-Monte Carlo of Lennard-Jonesium
 * Re-initialize checkpoint file
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

  FILE *fpi,*fpo;
  char fname[BSIZE];
  unsigned short seed[3];
  int n,nt,ntjob,ntprint,ntskip;
  int i;
  double disp,dr,rho,t;
  double fact,rhon;
  double *x,*y,*z;

  /* Read old checkpoint file */

  strcpy(fname,"mclj_out.dat");
  printf("         infile=[mclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  fread(&n,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  fread(&disp,sizeof(double),1,fpi);
  fread(&dr,sizeof(double),1,fpi);
  fread(&rho,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fclose(fpi);

  rhon=rho;

  /* User input */

  printf("              n=%13d\n",n);
  printf("            rho=[%12.5lf] ",rho);
  getval_d(&rhon);
  printf("              t=[%12.5lf] ",t);
  getval_d(&t);
  printf("           disp=[%12.5lf] ",disp);
  getval_d(&disp);
  printf("             dr=[%12.5lf] ",dr);
  getval_d(&dr);
  printf("         ntskip=[%12d] ",ntskip);
  getval_i(&ntskip);
  printf(" ntprint/ntskip=[%12d] ",ntprint);
  getval_i(&ntprint);
  printf("   ntjob/ntskip=[%12d] ",ntjob);
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

  /* Write new startup file */

  strcpy(fname,"mclj_in.dat");
  printf("        outfile=[ mclj_in.dat] ");
  getval_s(fname);

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
