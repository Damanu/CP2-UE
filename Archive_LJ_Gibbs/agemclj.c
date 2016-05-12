/**********************************************************************
 *
 * File: agemclj.c
 *
 * Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
 * Analyze checkpoint file
 *
 * 03-May-2010 (MN)
 * 19-APr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "getval.c"

#define max(a,b) ((a)>(b)?(a):(b))

#define NL '\n'
#define BSIZE 80

/**********************************************************************/

int main() {

/**********************************************************************/

  const double blue=0.0,green=0.0,radius=0.5,red=0.0;

  FILE *fpi,*fpo;
  char copy;
  char fname[BSIZE];
  unsigned short seed[3];
  int n,n1,nt,ntjob,ntprint,ntskip;
  int i,n2;
  double disp1,disp2,dv,ptr,pvm,t,v,v1;
  double c,c2;
  double *x,*y,*z;
  long long accrp1,accrp2,accrt,accrv,atryp1,atryp2,atryt1,atryt2,atryv,\
    an1,an2;
  double amu1,ap1,arho1,au1,av1,amu2,ap2,arho2,au2,av2;

  /* Read checkpoint file */

  strcpy(fname,"gemclj_out.dat");
  printf("         fname=[gemclj_out.dat] ");
  getval_s(fname);

  fpi=fopen(fname,"r");

  fread(&n,sizeof(int),1,fpi);
  fread(&n1,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  /* Check for zero-length run */

  if(nt<=0) {
    fclose(fpi);
    printf(" agemclj: empty file\n");
    exit(1);
  }

  /* Simulation parameters */

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

  /* Positions & accumulated averages */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(&accrp1,sizeof(long long),1,fpi);
  fread(&accrp2,sizeof(long long),1,fpi);
  fread(&accrt,sizeof(long long),1,fpi);
  fread(&accrv,sizeof(long long),1,fpi);
  fread(&atryp1,sizeof(long long),1,fpi);
  fread(&atryp2,sizeof(long long),1,fpi);
  fread(&atryt1,sizeof(long long),1,fpi);
  fread(&atryt2,sizeof(long long),1,fpi);
  fread(&atryv,sizeof(long long),1,fpi);
  fread(&amu1,sizeof(double),1,fpi);
  fread(&an1,sizeof(long long),1,fpi);
  fread(&ap1,sizeof(double),1,fpi);
  fread(&arho1,sizeof(double),1,fpi);
  fread(&au1,sizeof(double),1,fpi);
  fread(&av1,sizeof(double),1,fpi);
  fread(&amu2,sizeof(double),1,fpi);
  fread(&an2,sizeof(long long),1,fpi);
  fread(&ap2,sizeof(double),1,fpi);
  fread(&arho2,sizeof(double),1,fpi);
  fread(&au2,sizeof(double),1,fpi);
  fread(&av2,sizeof(double),1,fpi);

  fclose(fpi);

  /* Print results: simulation parameters */

  n2=n-n1;

  printf("\n");
  printf("            n1=%10d\n",n1);
  printf("            n2=%10d\n",n2);
  if(v<9999.9) {
    printf("             v=%10.5lf\n",v);
  } else {
    printf("             v=%12.5le\n",v);
  }
  printf("             t=%10.5lf\n",t);
  printf("         disp1=%10.5lf\n",disp1);
  printf("         disp2=%10.5lf\n",disp2);
  printf("            dv=%10.5lf\n",dv);
  printf("   prob(vmove)=%10.5lf\n",pvm);
  printf("   prob(trans)=%10.5lf\n",ptr);
  printf("\n");

  /* Averages */

  printf("            nt=%12d (*%5d)\n",nt,ntskip);
  printf("        accrp1=%12.5le\n",(double)accrp1/max(atryp1,1));
  printf("        accrp2=%12.5le\n",(double)accrp2/max(atryp2,1));
  printf("         accrv=%12.5le\n",(double)accrv/max(atryv,1));
  printf("         accrt=%12.5le\n",(double)accrt/max(atryt1+atryt2,1));
  printf("          <N1>=%12.5le\n",(double)an1/nt);
  printf("          <N2>=%12.5le\n",(double)an2/nt);
  printf("          <V1>=%12.5le\n",av1/nt);
  printf("          <V2>=%12.5le\n",av2/nt);
  printf("        <rho1>=%12.5le\n",arho1/nt);
  printf("        <rho2>=%12.5le\n",arho2/nt);
  printf("       <U1/N1>=%12.5le\n",au1/nt);
  printf("       <U2/N2>=%12.5le\n",au2/nt);
  printf("          <p1>=%12.5le\n",ap1/nt);
  printf("          <p2>=%12.5le\n",ap2/nt);
  if((amu1!=0.0)&&(atryt1!=0)) {
    printf("         <mu1>=%12.5le\n",-t*log(amu1/atryt1));
  }
  if((amu2!=0.0)&&(atryt2!=0)) {
    printf("         <mu2>=%12.5le\n",-t*log(amu2/atryt2));
  }
  printf("\n");

  /* Write PDB files? */

  printf(" Write PDB format to 'agemclj?.pdb'? [y] ");
  copy=(char)fgetc(stdin);

  if(copy=='y'||copy=='Y'||copy==NL) {

    /* Box 1 */

    c=cbrt(v1);
    c2=0.5*c;

    fpo=fopen("agemclj1.pdb","w");

    fprintf(fpo,"CRYST1"" %8.3lf %8.3lf %8.3lf",c,c,c);
    fprintf(fpo," %6.2f %6.2f %6.2f",90.0,90.0,90.0);
    fprintf(fpo,"  P1           1\n");

    for(i=0;i<=n1-1;i++) {
      fprintf(fpo,"HETATM%5d                    %7.3lf %7.3lf %7.3lf\n",\
              i+1,x[i]+c2,y[i]+c2,z[i]+c2);
    }
    
    fprintf(fpo,"COLOR ##### ####              ");
    fprintf(fpo," %7.3lf %7.3lf %7.3lf %5.2lf\n",red,green,blue,radius);
    fprintf(fpo,"END\n");


    fclose(fpo);

    /* Box 2 */

    c=cbrt(v-v1);
    c2=0.5*c;

    fpo=fopen("agemclj2.pdb","w");

    fprintf(fpo,"CRYST1"" %8.3lf %8.3lf %8.3lf",c,c,c);
    fprintf(fpo," %6.2f %6.2f %6.2f",90.0,90.0,90.0);
    fprintf(fpo,"  P1           1\n");

    for(i=n1;i<=n-1;i++) {
      fprintf(fpo,"HETATM%5d                    %7.3lf %7.3lf %7.3lf\n",\
              i+1,x[i]+c2,y[i]+c2,z[i]+c2);
    }
    
    fprintf(fpo,"COLOR ##### ####              ");
    fprintf(fpo," %7.3lf %7.3lf %7.3lf %5.2lf\n",red,green,blue,radius);
    fprintf(fpo,"END\n");

    fclose(fpo);

  }

  /* Deallocate arrays */

  free(x);
  free(y);
  free(z);

  return 0;

}

/**********************************************************************/
