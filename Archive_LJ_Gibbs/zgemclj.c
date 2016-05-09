/**********************************************************************
 *
 * File: zgemclj.c
 *
 * Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
 * Re-initialize checkpoint file
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

  FILE *fpi,*fpo;
  char fname[BSIZE];
  unsigned short seed[3];
  int n,n1,nt,ntjob,ntprint,ntskip;
  int i;
  double disp1,disp2,dv,ptr,pvm,t,v,v1;
  double fact,ptrold,pvmold,vold;
  double *x,*y,*z;

  /* Read old checkpoint file */

  strcpy(fname,"gemclj_out.dat");
  printf("         infile=[gemclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  fread(&n,sizeof(int),1,fpi);
  fread(&n1,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  fread(&disp1,sizeof(double),1,fpi);
  fread(&disp2,sizeof(double),1,fpi);
  fread(&dv,sizeof(double),1,fpi);
  fread(&ptr,sizeof(double),1,fpi);
  fread(&pvm,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);
  fread(&v,sizeof(double),1,fpi);
  fread(&v1,sizeof(double),1,fpi);

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

  ptrold=ptr;
  pvmold=pvm;
  vold=v;

  printf("              n=%15d\n",n);
  printf("             n1=%15d\n",n1);
  if(v<=99999999.9) {
    printf("              v=[%14.5lf] ",v);
  } else {
    printf("              v=[%14.5le] ",v);
  }
  getval_d(&v);
  printf("              t=[%14.5lf] ",t);
  getval_d(&t);
  printf("          disp1=[%14.5lf] ",disp1);
  getval_d(&disp1);
  printf("          disp2=[%14.5lf] ",disp2);
  getval_d(&disp2);
  printf("             dv=[%14.5lf] ",dv);
  getval_d(&dv);

  for(;;) {
    pvm=pvmold;
    ptr=ptrold;
    printf("    prob(vmove)=[%14.5lf] ",pvm);
    getval_d(&pvm);
    printf("    prob(trans)=[%14.5lf] ",ptr);
    getval_d(&ptr);
    if(ptr+pvm<1.0) {
      break;
    }
  }

  printf("         ntskip=[%14d] ",ntskip);
  getval_i(&ntskip);
  printf(" ntprint/ntskip=[%14d] ",ntprint);
  getval_i(&ntprint);
  printf("   ntjob/ntskip=[%14d] ",ntjob);
  getval_i(&ntjob);

  strcpy(fname,"gemclj_in.dat");
  printf("        outfile=[ gemclj_in.dat] ");
  getval_s(fname);

  /* Rescale positions */

  if(v!=vold) {
    v1=(v/vold)*v1;
    fact=cbrt(v/vold);
    for(i=0;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*x[i];
      z[i]=fact*x[i];
    }
  }

  /* Write new startup file */

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

