/**********************************************************************
 *
 * File: amclj.c
 *
 * NVT-Monte Carlo of Lennard-Jonesium
 * Analyze checkpoint file
 *
 * 12-May-2007 (MN)
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
  double disp,dr,rho,t;
  double c,c2,cv,fact,p,r,r1,r2;
  double *x,*y,*z;
  long long *ag;
  double accr,au,au2,aw;

  /* Read checkpoint file */

  strcpy(fname,"mclj_out.dat");
  printf(" fname=[mclj_out.dat] ");
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
    printf(" amclj: empty file\n");
    exit(0);
  }

  /* Simulation parameters */

  fread(&disp,sizeof(double),1,fpi);
  fread(&dr,sizeof(double),1,fpi);
  fread(&rho,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);

  /* Allocate arrays */

  c=cbrt(n/rho);
  c2=0.5*c;
  ndr=c2/dr;

  ag=(long long*)malloc(ndr*sizeof(long long));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions, accumulated averages & g(r) */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(&accr,sizeof(double),1,fpi);
  fread(&au,sizeof(double),1,fpi);
  fread(&au2,sizeof(double),1,fpi);
  fread(&aw,sizeof(double),1,fpi);

  fread(ag,sizeof(long long),ndr,fpi);

  fclose(fpi);

  /* Print results: simulation parameters */

  printf("\n");
  printf("     n=%9d\n",n);
  printf("   rho=%9.5lf\n",rho);
  printf("     t=%9.5lf\n",t);
  printf("  disp=%9.5lf\n",disp);
  printf("\n");

  /* Averages */

  cv=1.5+n*(au2/nt-(au/nt)*(au/nt))/(t*t);
  p=rho*(t-aw/(3*nt));
  printf("    nt=%12d (*%5d)\n",nt,ntskip);
  printf("  accr=%12.5e\n",accr/nt);
  printf(" <U>/N=%12.5e\n",au/nt);
  printf("  Cv/N=%12.5e\n",cv);
  printf("     p=%12.5e\n",p);
  printf("\n");

  /* Write g(r)? */

  printf(" Write g(r) to 'amclj.dat? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fact=3.0/(2.0*M_PI*rho*n*nt);
    fpo=fopen("amclj.dat","w");

    for(i=0;i<=ndr-2;i++) {
      r1=i*dr;
      r2=(i+1)*dr;
      r=0.5*(r1+r2);
      fprintf(fpo," %8.5lf %12.5e\n",\
	      r,fact*ag[i]/(r2*r2*r2-r1*r1*r1));
    }

    fclose(fpo);

  }

  if(copy!=NL) {
    copy=(char)fgetc(stdin);
  }
  
  /* Write PDB file? */

  printf(" Write PDB format to amclj.pdb? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fpo=fopen("amclj.pdb","w");

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
