/**********************************************************************
 *
 * File: agcmclj.c
 *
 * Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
 * Analyze checkpoint file
 *
 * 09-Apr-2010 (MN)
 * 19-Apr-2012
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
  double disp,dr,pcre,t,v,zz;
  double c,c2,chit,fact,p,r,r1,r2,rho,xmux;
  double *x,*y,*z;
  long long *ag;
  double accrc,accrd,accrm,an,an2,au,aw;

  /* Read checkpoint file */

  strcpy(fname,"gcmclj_out.dat");
  printf("      fname=[gcmclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname, "r");

  fread(&n,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  /* Check for zero-length run */

  if(nt<=0) {
    fclose(fpi);
    printf(" agcmclj: empty file\n");
    exit(1);
  }

  /* Simulation parameters */

  fread(&disp,sizeof(double),1,fpi);
  fread(&dr,sizeof(double),1,fpi);
  fread(&pcre,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);
  fread(&v,sizeof(double),1,fpi);
  fread(&zz,sizeof(double),1,fpi);

  /* Allocate arrays */

  ndr=0.5*cbrt(v)/dr;

  ag=(long long*)malloc(ndr*sizeof(long long));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions & accumulated averages */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(&accrc,sizeof(double),1,fpi);
  fread(&accrd,sizeof(double),1,fpi);
  fread(&accrm,sizeof(double),1,fpi);
  fread(&an,sizeof(double),1,fpi);
  fread(&an2,sizeof(double),1,fpi);
  fread(&au,sizeof(double),1,fpi);
  fread(&aw,sizeof(double),1,fpi);

  fread(ag,sizeof(long long),ndr,fpi);

  fclose(fpi);

  /* Print results: simulation parameters */

  printf("\n");
  printf("          t=%14.7lf\n",t);
  printf("          v=%14.7lf\n",v);
  printf("          z=%14.7lf\n",zz);
  printf("       disp=%14.7lf\n",disp);
  printf(" p(cre/del)=%14.7lf\n",pcre);
  printf("\n");

  /* Averages */

  rho=an/(v*nt);
  xmux=t*log(zz/rho);
  p=(an*t-aw/3.0)/(v*nt);
  chit=(an2/an-an/nt)/(rho*t);

  printf("         nt=%14d (*%5d)\n",nt,ntskip);
  printf("      accrc=%14.7le\n",accrc/nt);
  printf("      accrd=%14.7le\n",accrd/nt);
  printf("      accrm=%14.7le\n",accrm/nt);
  printf("      mu_ex=%14.7le\n",xmux);
  printf("        <n>=%14.7le\n",an/nt);
  printf("      <rho>=%14.7le\n",rho);
  printf("    <u>/<n>=%14.7le\n",au/an);
  printf("          p=%14.7le\n",p);
  printf("      chi_t=%14.7le\n",chit);
  printf("\n");

  /* Write g(r) to file? */

  printf("       Write g(r) to 'agcmclj.dat'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fact=3.0/(2.0*M_PI*rho*an);

    fpo=fopen("agcmclj.dat","w");

    for(i=0;i<=ndr-1;i++) {
      r1=i*dr;
      r2=(i+1)*dr;
      r=0.5*(r1+r2);
      fprintf(fpo," %8.5f %14.7le\n",r,fact*ag[i]/(r2*r2*r2-r1*r1*r1));
    }

    fclose(fpo);

  }

  if(copy!=NL) {
    copy=(char)fgetc(stdin);
  }

  /* Write PDB file? */

  printf(" Write PDB format to 'agcmclj.pdb'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    c=cbrt(v);
    c2=0.5*c;

    fpo=fopen("agcmclj.pdb","w");

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
  free(x);
  free(y);
  free(z);

  return 0;

}

/**********************************************************************/
