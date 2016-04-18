/**********************************************************************
 *
 * File: nptmclj.c
 *
 * NpT-Monte Carlo of Lennard-Jonesium
 *
 * 03-Apr-2010 (MN)
 * 19-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

/* Parameters & global variables */

double alnmax;
unsigned short seed[3];
int n,nt,ntjob,ntprint,ntskip;
int naccp,naccv,ndr,ntryp,ntryv;
double disp,dr,dlnv,p,pvm,t,v;
double c,c2,c2m1,drm1,r2max,rc,rc2,su0;
double *x,*y,*z;
double *uattb,*urepb,*uattn,*urepn;
double **uatt,**urep;
long long *ag;
double accrp,accrv,arho,au,aupv2,av,av2;

/**********************************************************************/

void getcf() {

/**********************************************************************/

  FILE *fpi;
  int i;

  /* Maximum argument for exponential */

  alnmax=log(0.1*DBL_MAX);

  /* Read startup/checkpoint file */

  fpi=fopen("nptmclj_in.dat","r");

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

  ndr=0.5*cbrt(v)/dr;

  ag=(long long*)malloc(ndr*sizeof(long long));
  uatt=(double**)malloc(n*sizeof(double*));
  urep=(double**)malloc(n*sizeof(double*));
  uattb=(double*)malloc(n*n*sizeof(double));
  urepb=(double*)malloc(n*n*sizeof(double));
  uattn=(double*)malloc(n*sizeof(double));
  urepn=(double*)malloc(n*sizeof(double));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  /* Read accumulated averages/clear accumulators */

  if(nt>0) {

    fread(&accrp,sizeof(double),1,fpi);
    fread(&accrv,sizeof(double),1,fpi);
    fread(&arho,sizeof(double),1,fpi);
    fread(&au,sizeof(double),1,fpi);
    fread(&aupv2,sizeof(double),1,fpi);
    fread(&av,sizeof(double),1,fpi);
    fread(&av2,sizeof(double),1,fpi);
    fread(&ndr,sizeof(int),1,fpi);
    fread(ag,sizeof(long long),ndr,fpi);

  } else {

    accrp=0.0;
    accrv=0.0;
    arho=0.0;
    au=0.0;
    aupv2=0.0;
    av=0.0;
    av2=0.0;

    for(i=0;i<=ndr-1;i++) {
      ag[i]=0;
    }
    
  }

  fclose(fpi);

  /* Box parameters */

  c=cbrt(v);
  c2=0.5*c;
  c2m1=1.0/c2;
  r2max=c2*c2;
  drm1=1.0/dr;

  /* Cutoff & tail correction */

  rc=c2;
  rc2=rc*rc;
  su0=2.0*M_PI*n*n/v*(4.0/(9.0*pow(rc,9))-4.0/(3.0*pow(rc,3)));

  /* Initialize RNG & acceptance counters */

  seed48(seed);

  ntryp=0;
  ntryv=0;
  naccp=0;
  naccv=0;

  return;

}

/**********************************************************************/

void uinit() {

/**********************************************************************/

  int i,j;
  double dx,dy,dz,r2,rm6;

  /* Initialize pair interaction matrix
   *
   * urep = 4/r^12, uatt = -4/r^6 */

  for(i=0;i<=n-1;i++) {
    urep[i]=urepb+i*n;
    uatt[i]=uattb+i*n;
  }

  for(i=0;i<=n-2;i++) {
    urep[i][i]=0.0;
    uatt[i][i]=0.0;
    for(j=i+1;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
	rm6=1.0/(r2*r2*r2);
	urep[i][j]=4.0*rm6*rm6;
	uatt[i][j]=-4.0*rm6;
      } else {
	urep[i][j]=0.0;
	uatt[i][j]=0.0;
      }
      urep[j][i]=urep[i][j];
      uatt[j][i]=uatt[i][j];
    }
  }
  urep[n-1][n-1]=0.0;
  uatt[n-1][n-1]=0.0;

  return;

}

/**********************************************************************/

void pmove() {

/**********************************************************************/

  int accept;
  int i,j;
  double dx,dy,dz,r2,rm6,su,xin,yin,zin;

  ntryp=ntryp+1;

  /* Random particle */

  i=min((int)(drand48()*n),n-1);

  /* Trial move */

  xin=x[i]+disp*(drand48()-0.5);
  xin=xin-(int)(xin*c2m1)*c;
  yin=y[i]+disp*(drand48()-0.5);
  yin=yin-(int)(yin*c2m1)*c;
  zin=z[i]+disp*(drand48()-0.5);
  zin=zin-(int)(zin*c2m1)*c;

  /* Energy difference */

  su=0.0;
  for(j=0;j<=n-1;j++) {
    if(j==i) {
      urepn[j]=0.0;
      uattn[j]=0.0;
    } else {
      dx=x[j]-xin;
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-yin;
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-zin;
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
	rm6=1.0/(r2*r2*r2);
	urepn[j]=4.0*rm6*rm6;
	uattn[j]=-4.0*rm6;
      } else {
	urepn[j]=0.0;
	uattn[j]=0.0;
      }
      su=su+((urepn[j]+uattn[j])-(urep[i][j]+uatt[i][j]));
    }
  }

  /* Acceptance test */

  if(su<=0.0) {
    accept=TRUE;
  } else {
    if(su/t<alnmax) {
      if(drand48()<=exp(-su/t)) {
	accept=TRUE;
      } else {
	accept=FALSE;
      }
    } else {
      accept=FALSE;
    }
  }

  /* Accept move & update pair interaction matrix */

  if(accept) {
    naccp=naccp+1;
    x[i]=xin;
    y[i]=yin;
    z[i]=zin;
    for(j=0;j<=n-1;j++) {
      urep[i][j]=urepn[j];
      urep[j][i]=urepn[j];
      uatt[i][j]=uattn[j];
      uatt[j][i]=uattn[j];
    }
  }

  return;

}

/**********************************************************************/

void vmove() {

/**********************************************************************/

  int accept;
  int i,j;
  double fact,fatt,frep,rcn,su,sun0,vn;

  ntryv=ntryv+1;

  /* Volume change */

  vn=v*exp(dlnv*(drand48()-0.5));

  /* New cutoff & tail correction */

  rcn=0.5*cbrt(vn);
  sun0=2.0*M_PI*n*n/vn*(4.0/(9.0*pow(rcn,9))-4.0/(3.0*pow(rcn,3)));

  /* "Energy" difference */

  frep=pow(v/vn,4)-1.0;
  fatt=pow(v/vn,2)-1.0;
  su=sun0-su0;
  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      su=su+(frep*urep[i][j]+fatt*uatt[i][j]);
    }
  }
  su=su+p*(vn-v)-(n+1)*t*log(vn/v);

  /* Acceptance test */

  if(su<=0.0) {
    accept=TRUE;
  } else {
    if(su/t<alnmax) {
      if(drand48()<=exp(-su/t)) {
	accept=TRUE;
      } else {
	accept=FALSE;
      }
    } else {
      accept=FALSE;
    }
  }

  /* Accept */

  if(accept) {
    naccv=naccv+1;

    /* Positions */

    fact=cbrt(vn/v);
    for(i=0;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }

    /* Box parameters */

    v=vn;
    c=cbrt(v);
    c2=0.5*c;
    c2m1=1.0/c2;
    r2max=c2*c2;
    ndr=min(ndr,(int)(c2/dr));

    /* Cutoff & tail correction */

    rc=rcn;
    rc2=rc*rc;
    su0=sun0;

    /* Update pair interaction matrix */

    uinit();

  }

  return;

}

/**********************************************************************/

void means() {

/**********************************************************************/

  int i,j,k;
  double dx,dy,dz,r2,su;

  /* Potential energy & g(r) */

  su=su0;
  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      su=su+(urep[i][j]+uatt[i][j]);
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<r2max) {
	k=(int)(sqrt(r2)*drm1);
	if(k<=ndr-1) {
	  ag[k]=ag[k]+1;
	}
      }
    }
  }
  su=su/n;

  /* Print control variables */

  if((ntprint>0)&&(nt%ntprint==0)) {
    printf(" %10d %12.5e %12.5e %12.5e %12.5e %12.5e\n",\
	   nt,(double)naccp/max(ntryp,1),(double)naccv/max(ntryv,1),\
	   su,v,n/v);
  }

  /* Accumulate averages */

  accrp=accrp+(double)naccp/max(ntryp,1);
  accrv=accrv+(double)naccv/max(ntryv,1);
  arho=arho+(n/v);
  au=au+su;
  aupv2=aupv2+pow(su+p*v/n,2);
  av=av+v;
  av2=av2+v*v;

  /* Clear acceptance counters */

  naccp=0;
  naccv=0;
  ntryp=0;
  ntryv=0;

  return;

}

/**********************************************************************/

void putcf() {

/**********************************************************************/

  FILE *fpo;
  unsigned short *lastseed;

  /* RNG seed */

  lastseed=seed48(seed);
  seed[0]=lastseed[0];
  seed[1]=lastseed[1];
  seed[2]=lastseed[2];

  /* Write checkpoint file */

  fpo=fopen("nptmclj_out.dat","w");

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

  fwrite(&accrp,sizeof(double),1,fpo);
  fwrite(&accrv,sizeof(double),1,fpo);
  fwrite(&arho,sizeof(double),1,fpo);
  fwrite(&au,sizeof(double),1,fpo);
  fwrite(&aupv2,sizeof(double),1,fpo);
  fwrite(&av,sizeof(double),1,fpo);
  fwrite(&av2,sizeof(double),1,fpo);

  fwrite(&ndr,sizeof(int),1,fpo);
  fwrite(ag,sizeof(long long),ndr,fpo);

  fclose(fpo);

  /* Deallocate arrays */

  free(ag);
  free(uatt);
  free(urep);
  free(uattb);
  free(urepb);
  free(uattn);
  free(urepn);
  free(x);
  free(y);
  free(z);

  return;

}


/**********************************************************************/

int main() {

/**********************************************************************/

  int i,j;

  /* Read startup/checkpoint file */

  getcf();

  /* Initialize pair interaction matrix */

  uinit();

  /* Do ntskip*ntjob passes/volume changes */

  for(i=1;i<=ntjob;i++) {
    nt=nt+1;
    for(j=1;j<=ntskip*n;j++) {
      if(drand48()<pvm) {
	vmove();            /* Volume change */
      } else {
	pmove();            /* Particle move */
      }
    }
    means();
  }

  /* Write checkpoint file */

  putcf();

  return 0;

}

/**********************************************************************/
