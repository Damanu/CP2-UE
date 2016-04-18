/**********************************************************************
 *
 * File: znptmclj.c
 *
 * NpT-Monte Carlo of Lennard-Jonesium
 * Re-initialize checkpoint file
 *
 * 03-Apr-2010 (MN)
 * 17-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  double disp,dlnv,dr,p,pvm,t,v;
  double *x,*y,*z;

  /* Read old checkpoint file */

  strcpy(fname,"nptmclj_out.dat");
  printf("         infile=[nptmclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  fread(&n,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  fread(&disp,sizeof(double),1,fpi);
  fread(&dlnv,sizeof(double),1,fpi);
  fread(&dr,sizeof(double),1,fpi);
  fread(&p,sizeof(double),1,fpi);
  fread(&pvm,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);
  fread(&v,sizeof(double),1,fpi);

  /* Allocate arrays */

  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fclose(fpi);

  /* User input */

  printf("              n=%16d\n",n);
  printf("              p=[%15.5lf] ",p);
  getval_d(&p);
  printf("              t=[%15.5lf] ",t);
  getval_d(&t);
  printf("           disp=[%15.5lf] ",disp);
  getval_d(&disp);
  printf("           dv/v=[%15.5lf] ",dlnv);
  getval_d(&dlnv);
  printf("    prob(vmove)=[%15.5lf] ",pvm);
  getval_d(&pvm);
  printf("             dr=[%15.5lf] ",dr);
  getval_d(&dr);
  printf("         ntskip=[%15d] ",ntskip);
  getval_i(&ntskip);
  printf(" ntprint/ntskip=[%15d] ",ntprint);
  getval_i(&ntprint);
  printf("          ntjob=[%15d] ",ntjob);
  getval_i(&ntjob);

  /* Write new startup file */

  strcpy(fname,"nptmclj_in.dat");
  printf("        outfile=[nptmclj_out.dat] ");
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
