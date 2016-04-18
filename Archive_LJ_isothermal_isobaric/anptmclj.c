/**********************************************************************
 *
 * File: anptmclj.c
 *
 * NpT-Monte Carlo of Lennard-Jonesium
 * Analyze checkpoint file
 *
 * 03-Apr-2010 (MN)
 * 23-Apr-2012
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

  const double blue=0.0,green=0.0,radius=0.5,red=1.0;

  FILE *fpi,*fpo;
  char copy;
  char fname[BSIZE];
  unsigned short seed[3];
  int n,nt,ntjob,ntprint,ntskip;
  int i,ndr;
  double disp,dlnv,dr,p,pvm,t,v;
  double c,c2,chit,cp,fact,r,r1,r2;
  double *x,*y,*z;
  long long *ag;
  double accrp,accrv,arho,au,aupv2,av,av2;

  /* Read checkpoint file */

  strcpy(fname,"nptmclj_out.dat");
  printf(" fname=[nptmclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  fread(&n,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  /* Check for zero-length run */

  if(nt<=0) {
    fclose(fpi);
    printf(" anptmclj: empty file\n");
    exit(1);
  }

  /* Simulation parameters */

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

  /* Positions & accumulated averages */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(&accrp,sizeof(double),1,fpi);
  fread(&accrv,sizeof(double),1,fpi);
  fread(&arho,sizeof(double),1,fpi);
  fread(&au,sizeof(double),1,fpi);
  fread(&aupv2,sizeof(double),1,fpi);
  fread(&av,sizeof(double),1,fpi);
  fread(&av2,sizeof(double),1,fpi);

  /* g(r) */

  fread(&ndr,sizeof(int),1,fpi);
  ag=(long long*)malloc(ndr*sizeof(long long));
  fread(ag,sizeof(long long),ndr,fpi);
  
  fclose(fpi);

  /* Print results: simulation parameters */

  printf("\n");
  printf("           n=%9d\n",n);
  printf("           p=%9.5lf\n",p);
  printf("           t=%9.5lf\n",t);
  printf("        disp=%9.5lf\n",disp);
  printf("        dv/v=%9.5lf\n",dlnv);
  printf(" prob(vmove)=%9.5lf\n",pvm);
  printf("\n");

  /* Averages */

  cp=1.5+n*(aupv2/nt-pow((au+p*av/n)/nt,2))/pow(t,2);
  chit=(av2/nt-pow(av/nt,2))/(av/nt)/t;

  printf("          nt=%12d (*%5d)\n",nt,ntskip);
  printf("       accrp=%12.5e\n",accrp/nt);
  printf("       accrv=%12.5e\n",accrv/nt);
  printf("       <rho>=%12.5e\n",arho/nt);
  printf("       <U>/N=%12.5e\n",au/nt);
  printf("       <V>/N=%12.5e\n",av/nt);
  printf("        Cp/N=%12.5e\n",cp);
  printf("       chi_t=%12.5e\n",chit);
  printf("\n");

  /* Write g(r) to file? */

  printf(" Write g(r) to 'amclj.dat? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fact=3.0/(2.0*M_PI*arho*n);

    fpo=fopen("anptmclj.dat","w");

    for(i=0;i<=ndr-1;i++) {
      r1=i*dr;
      r2=(i+1)*dr;
      r=0.5*(r1+r2);
      fprintf(fpo," %8.5f %12.5e\n",r,fact*ag[i]/(r2*r2*r2-r1*r1*r1));
    }

    fclose(fpo);

  }

  if(copy!=NL) {
    copy=(char)fgetc(stdin);
  }

  /* Write PDB file?*/

  printf(" Write PDB format to 'anptmclj.pdb'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    c=cbrt(v);
    c2=0.5*c;

    fpo=fopen("anptmclj.pdb","w");

    fprintf(fpo,"CRYST1"" %8.3lf %8.3lf %8.3lf",c,c,c);
    fprintf(fpo," %6.2f %6.2f %6.2f",90.0,90.0,90.0);
    fprintf(fpo,"  P1           1\n");

    for(i=0;i<=n-1;i++) {
      fprintf(fpo,"HETATM%5d                    %7.3lf %7.3lf %7.3lf\n",\
              i+1,x[i]+c2,y[i]+c2,z[i]+c2);
    }

    fprintf(fpo,"COLOR ##### ####              ");
    fprintf(fpo," %7.3lf %7.3lf %7.3lf %5.2lf\n",red,green,blue,radius);
    fprintf(fpo,"END\n");

    fclose(fpo);

  }

  /* Deallocate arrays */

  free(ag);
/*  free(x);
  free(y);
  free(z); */

  return 0;

}

/**********************************************************************/
