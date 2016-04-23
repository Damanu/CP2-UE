/***********************************************************************
 *
 * File: zgcmclj.c
 *
 * Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
 * Re-initialize checkpoint file
 *
 * 08-Apr-2010 (MN)
 * 18-Apr-2012
 *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.h"

#define BSIZE 80

/***********************************************************************/

int main() {

/***********************************************************************/

  FILE *fpi,*fpo;
  char fname[BSIZE];
  unsigned short seed[3];
  int n,nt,ntjob,ntprint,ntskip;
  int i;
  double disp,dr,pcre,t,v,zz;
  double fact,v_old;
  double *x,*y,*z;

  /* Read old checkpoint file */

  strcpy(fname,"gcmclj_out.dat");
  printf("         infile=[gcmclj_out.dat] ");
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
  fread(&pcre,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);
  fread(&v,sizeof(double),1,fpi);
  fread(&zz,sizeof(double),1,fpi);

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

  printf("              t=[%14.7lf] ",t);
  getval_d(&t);
  v_old=v;
  printf("              v=[%14.7lf] ",v);
  getval_d(&v);
  printf("              n=%15d\n",n);
  printf("              z=[%14.7lf] ",zz);
  getval_d(&zz);

  for(;;) {
    printf("     p(cre/del)=[%14.7lf] ",pcre);
    getval_d(&pcre);
    if(pcre<0.5) {
      break;
    }
  }

  printf("           disp=[%14.7lf] ",disp);
  getval_d(&disp);
  printf("             dr=[%14.7lf] ",dr);
  getval_d(&dr);
  printf("         ntskip=[%14d] ",ntskip);
  getval_i(&ntskip);
  printf(" ntprint/ntskip=[%14d] ",ntprint);
  getval_i(&ntprint);
  printf("   ntjob/ntskip=[%14d] ",ntjob);
  getval_i(&ntjob);

  strcpy(fname,"gcmclj_in.dat");
  printf("        outfile=[ gcmclj_in.dat] ");
  getval_s(fname);

  /* Rescale positions */

  if(v!=v_old) {
    fact=cbrt(v/v_old);
    for(i=0;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }
  }

  /* Write new startup file */

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
  fwrite(&pcre,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);
  fwrite(&v,sizeof(double),1,fpo);
  fwrite(&zz,sizeof(double),1,fpo);

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

/***********************************************************************/
