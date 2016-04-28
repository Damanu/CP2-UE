/**********************************************************************
 *
 * File: aisomdlj.c
 *
 * Gaussian isokinetic (NVT) Molecular Dynamics of Lennard-Jonesium
 * Analyze checkpoint file
 *
 * 20-Nov-2007 (MN)
 * 05-May-2012
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
  int i,ndr;
  int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
  long long *ag;
  double c,c2,cv,fact,p,r,r1,r2;
  double dr,dt,rho,t;
  double *x,*y,*z;
  double ak,au,au2,aw;
  double *acf;

  /* Read checkpoint file */

  strcpy(fname,"isomdlj_out.dat");
  printf(" fname=[isomdlj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  /* Simulation parameters */

  fread(&n,sizeof(int),1,fpi);
  fread(&ncor,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntaskip,sizeof(int),1,fpi);
  fread(&ntcskip,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntorig,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);

  fread(&dr,sizeof(double),1,fpi);
  fread(&dt,sizeof(double),1,fpi);
  fread(&rho,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);

  /* Check for zero-length run */

  if(nt<=0) {
    fclose(fpi);
    printf(" aisomdlj: empty file\n");
    exit(1);
  }

  /* Allocate arrays */

  ndr=(int)(0.5*cbrt(n/rho)/dr);
  
  acf=(double*)malloc((ncor+1)*sizeof(double));
  ag=(long long*)malloc(ndr*sizeof(long long));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  /* Skip velocities: vx,vy,vz */

  fseek(fpi,(long)(3*n*sizeof(double)),SEEK_CUR);

  /* Accumulated averages */

  fread(&ak,sizeof(double),1,fpi);
  fread(&au,sizeof(double),1,fpi);
  fread(&au2,sizeof(double),1,fpi);
  fread(&aw,sizeof(double),1,fpi);

  /* Pair distribution function g(r) */

  fread(ag,sizeof(long long),ndr,fpi);

  /* Velocity autocorrelation function C(t) */

  if(nt/ntcskip>ncor) {
    /* Skip ring buffer */
    for(i=0;i<=ncor;i++) {
      fseek(fpi,(long)(3*n*sizeof(double)),SEEK_CUR); /* vxt,vyt,vzt */
    }
    fread(acf,sizeof(double),ncor+1,fpi);             /* C(t) */
  }

  fclose(fpi);

  /* Print results: simulation parameters */

  printf("\n");
  printf("              n=%12d\n",n);
  printf("            rho=%12.5lf\n",rho);
  printf("              t=%12.5lf\n",t);
  printf("             dt=%12.5lf\n",dt);
  printf("        ntaskip=%12d\n",ntaskip);
  printf("        ntcskip=%12d\n",ntcskip);
  printf("   ncor/ntcskip=%12d\n",ncor);
  printf(" ntorig/ntcskip=%12d\n",ntorig);

  /* Averages */

  cv=1.5+n*(au2/(nt/ntaskip)-pow(au/(nt/ntaskip),2.0))/\
    pow(ak/(1.5*(nt/ntaskip)),2.0);
  p=rho*(2.0*ak-aw)/(3.0*(nt/ntaskip));

  printf("\n");
  printf("    nt=%12d\n",nt);
  printf("   <t>=%12.5le\n",ak/(1.5*(nt/ntaskip)));
  printf(" <k>/n=%12.5le\n",ak/(nt/ntaskip));
  printf(" <u>/n=%12.5le\n",au/(nt/ntaskip));
  printf("     p=%12.5le\n",p);
  printf("  cv/n=%12.5le\n",cv);
  printf("\n");

  /* Write g(r) to file? */

  printf(" Write g(r) to 'aisomdlj1.dat'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fpo=fopen("aisomdlj1.dat","w");

    fact=3.0/(2.0*M_PI*rho*n*(nt/ntaskip));
    for(i=0;i<=ndr-1;i++) {
      r1=i*dr;
      r2=(i+1)*dr;
      r=0.5*(r1+r2);
      fprintf(fpo," %8.5lf %12.5le\n",\
	      r,fact*ag[i]/(r2*r2*r2-r1*r1*r1));
    }

    fclose(fpo);

  }

  if(copy!=NL) {
    copy=(char)fgetc(stdin);
  }

  /* Write C(t) to file? */

  if(nt/ntcskip>ncor) {

    printf(" Write c(t) to 'aisomdlj2.dat'? [y] ");
    copy=(char)fgetc(stdin);

    if(copy=='y'||copy=='Y'||copy==NL) {

      fpo=fopen("aisomdlj2.dat","w");

      for(i=0;i<=ncor;i++) {
	fprintf(fpo," %8.5lf %12.5e %12.5e\n",i*ntcskip*dt,\
		acf[i]/((nt/ntcskip-(ncor+1))/ntorig+1), /* Unnormalized */\
		acf[i]/acf[0]);                          /* Normalized   */
      }

      fclose(fpo);

    }

  }

  if(copy!=NL) {
    copy=(char)fgetc(stdin);
  }

  /* Write PDB file? */

  printf(" Write PDB format to 'aisomdlj.pdb'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    fpo=fopen("aisomdlj.pdb","w");

    c=cbrt(n/rho);
    c2=0.5*c;
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

  free(acf);
  free(ag);
  free(x);
  free(y);
  free(z);

  return 0;

}

/**********************************************************************/
