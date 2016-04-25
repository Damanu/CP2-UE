/**********************************************************************
 *
 * File: gcmclj.c
 *
 * Grand Canonical (TVmu) Monte Carlo of Lennard-Jonesium
 *
 * 09-Apr-2010 (MN)
 * 19-Apr-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define TRUE 1
#define FALSE 0
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

/* Parameters & global variables */

int n,nt,ntjob,ntprint,ntskip;
int naccc,naccd,naccm,ndr,nmax;
double disp,dr,pcre,t,v,zz;
double alnmax,c,c2,c2m1,drm1,r2max,rc,rc2,su0,sw0,zv;
double *x,*y,*z;
long long *ag;
double accrc,accrd,accrm,an,an2,au,aw;

/**********************************************************************/

void getcf() {

/**********************************************************************/

  FILE *fpi;
  unsigned short seed[3];
  int i;

  /* Maximum argument for exponential */

  alnmax=log(0.1*DBL_MAX);

  /* Read startup/checkpoint file */

  fpi=fopen("gcmclj_in.dat","r");

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

  ndr=0.5*cbrt(v)/dr;

  ag=(long long*)malloc(ndr*sizeof(long long));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  /* Read accumulated averages/clear accumulators */

  if(nt>0) {

    fread(&accrc,sizeof(double),1,fpi);
    fread(&accrd,sizeof(double),1,fpi);
    fread(&accrm,sizeof(double),1,fpi);
    fread(&an,sizeof(double),1,fpi);
    fread(&an2,sizeof(double),1,fpi);
    fread(&au,sizeof(double),1,fpi);
    fread(&aw,sizeof(double),1,fpi);

    fread(ag,sizeof(long long),ndr,fpi);

  } else {

    accrc=0.0;
    accrd=0.0;
    accrm=0.0;
    an=0.0;
    an2=0.0;
    au=0.0;
    aw=0.0;

    for(i=0;i<=ndr-1;i++) {
      ag[i]=0;
    }

  }

  fclose(fpi);

  /* Parameters */

  nmax=n;

  c=cbrt(v);
  c2=0.5*c;
  c2m1=1.0/c2;
  drm1=1.0/dr;
  r2max=c2*c2;
  zv=zz*v;

  /* Cutoff & tail corrections */

  rc=0.5*c;
  rc2=rc*rc;

  su0=2.0*M_PI/v*(4.0/(9.0*pow(rc,9))-4.0/(3.0*pow(rc,3)));
  sw0=2.0*M_PI/v*(24.0/(3.0*pow(rc,3))-48.0/(9.0*pow(rc,9)));

  /* Initialize RNG */

  seed48(seed);

  /* Clear acceptance counters */

  naccc=0;
  naccd=0;
  naccm=0;

  return ;

}

/**********************************************************************/

void create() {

/**********************************************************************/

  int accept;
  int i,j;
  double dx,dy,dz,r2,rm6,su,xn,yn,zn;
  double *xt,*yt,*zt;

  /* Random position */

  xn=(drand48()-0.5)*c;
  yn=(drand48()-0.5)*c;
  zn=(drand48()-0.5)*c;

  /* Energy difference */

  su=(2*n+1)*su0;

  for(j=0;j<=n-1;j++) {
    dx=x[j]-xn;
    dx=dx-(int)(dx*c2m1)*c;
    dy=y[j]-yn;
    dy=dy-(int)(dy*c2m1)*c;
    dz=z[j]-zn;
    dz=dz-(int)(dz*c2m1)*c;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<rc2) {
      rm6=1.0/(r2*r2*r2);
      su=su+(4.0*rm6-4.0)*rm6;
    }
  }

  /* Acceptance test */

  if(su/t>alnmax) {
    accept=FALSE;
  } else {
    if(zv*exp(-su/t)/(n+1)>=1.0) {
      accept=TRUE;
    } else {
      if(drand48()<=zv*exp(-su/t)/(n+1)) {
	accept=TRUE;
      } else {
	accept=FALSE;
      }
    }
  }

  /* Accept move & append particle to arrays */

  if(accept) {

    naccc=naccc+1;

    /* Extend arrays */

    if(n>=nmax) {
      xt=(double*)malloc(n*sizeof(double));
      yt=(double*)malloc(n*sizeof(double));
      zt=(double*)malloc(n*sizeof(double));

      for(i=0;i<=n-1;i++) {
	xt[i]=x[i];
	yt[i]=y[i];
	zt[i]=z[i];
      }
      free(x);
      free(y);
      free(z);

      nmax=n+1;
      x=(double*)malloc(nmax*sizeof(double));
      y=(double*)malloc(nmax*sizeof(double));
      z=(double*)malloc(nmax*sizeof(double));

      for(i=0;i<=n-1;i++) {
	x[i]=xt[i];
	y[i]=yt[i];
	z[i]=zt[i];
      }
      free(xt);
      free(yt);
      free(zt);

    }

    x[n]=xn;
    y[n]=yn;
    z[n]=zn;
    n=n+1;

  }

  return;

}

/**********************************************************************/

void delete() {

/**********************************************************************/

  int accept;
  int i,j;
  double dx,dy,dz,r2,rm6,su;

  /* Empty box? */

  if(n<1) {
    return;
  }

  /* Random particle */

  i=min((int)(drand48()*n),n-1);

  /* Energy difference (ok for n=1) */

  su=-(2*n-1)*su0;

  for(j=0;j<=n-1;j++) {
    if(j!=i) {
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
	rm6=1.0/(r2*r2*r2);
	su=su-(4.0*rm6-4.0)*rm6;
      }
    }
  }

  /* Acceptance test */

  if(su/t>alnmax) {
    accept=FALSE;
  } else {
    if(n*exp(-su/t)/zv>=1.0) {
      accept=TRUE;
    } else {
      if(drand48()<=n*exp(-su/t)/zv) {
	accept=TRUE;
      } else {
	accept=FALSE;
      }
    }
  }

  /* Accept move & eliminate particle from arrays */

  if(accept) {
    naccd=naccd+1;
    x[i]=x[n-1];
    y[i]=y[n-1];
    z[i]=z[n-1];
    n=n-1;
  }

  return;

}

/**********************************************************************/

void move() {

/**********************************************************************/

  int accept;
  int i,j;
  double dx,dy,dz,r2,rm6,su,xn,yn,zn;

  /* Select random particle */

  i=min((int)(drand48()*n),n-1);

  /* Trial move (displacement) */

  xn=x[i]+(drand48()-0.5)*disp;
  xn=xn-(int)(xn*c2m1)*c;
  yn=y[i]+(drand48()-0.5)*disp;
  yn=yn-(int)(yn*c2m1)*c;
  zn=z[i]+(drand48()-0.5)*disp;
  zn=zn-(int)(zn*c2m1)*c;

  /* Energy difference */

  su=0.0;

  for(j=0;j<=n-1;j++) {
    if(j!=i) {

      /* New position */

      dx=x[j]-xn;
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-yn;
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-zn;
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
	rm6=1.0/(r2*r2*r2);
	su=su+(4.0*rm6-4.0)*rm6;
      }

      /* Old position */

      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
	rm6=1.0/(r2*r2*r2);
	su=su-(4.0*rm6-4.0)*rm6;
      }

    }
  }

  /* Acceptance test */

  if(su<=0.0) {
    accept=TRUE;
  } else {
    if(su/t>alnmax) {
      accept=FALSE;
    } else {
      if(drand48()<=exp(-su/t)) {
	accept=TRUE;
      } else {
	accept=FALSE;
      }
    }
  }

  /* Accept move */

  if(accept) {
    naccm=naccm+1;
    x[i]=xn;
    y[i]=yn;
    z[i]=zn;
  }

  return;

}

/**********************************************************************/

void means() {

/**********************************************************************/

  int i,ir,j;
  double dx,dy,dz,r2,rm6,su,sw,p;

  /* Potential energy, virial & g(r) */

  su=n*n*su0;
  sw=n*n*sw0;

  for(i=0;i<=n-2;i++) {
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
	su=su+(4.0*rm6-4.0)*rm6;
	sw=sw+(24.0-48.0*rm6)*rm6;
      }

      if(r2<r2max) {
	ir=(int)(drm1*sqrt(r2));
	if(ir<ndr) {
	  ag[ir]=ag[ir]+1;
	}
      }

    }
  }

  /* Print control variables */

  if((ntprint>0)&&(nt%ntprint==0)) {
    p=n/v*(t-sw/3.0);
    printf("%d %lf %lf %lf %d %lg %lg %lg\n",nt,\
	   naccc/max(pcre*n*ntskip,1.0),\
	   naccd/max(pcre*n*ntskip,1.0),\
	   naccm/max((1.0-2.0*pcre)*n*ntskip,1.0),\
	   n,n/v,su/n,p);
  }

  /* Accumulate averages */

  accrc=accrc+naccc/max(pcre*n*ntskip,1.0);
  accrd=accrd+naccd/max(pcre*n*ntskip,1.0);
  accrm=accrm+naccm/max((1.0-2.0*pcre)*n*ntskip,1.0);
  an=an+n;
  an2=an2+n*n;
  au=au+su;
  aw=aw+sw;

  /* Clear acceptance counters */

  naccc=0;
  naccd=0;
  naccm=0;

  return;

}

/**********************************************************************/

void putcf() {

/**********************************************************************/

  FILE *fpo;
  unsigned short dummy[3];
  unsigned short *seed;

  /* RNG seed */

  seed=seed48(dummy);

  /* Write checkpoint file */

  fpo=fopen("gcmclj_out.dat","w");

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

  fwrite(&accrc,sizeof(double),1,fpo);
  fwrite(&accrd,sizeof(double),1,fpo);
  fwrite(&accrm,sizeof(double),1,fpo);
  fwrite(&an,sizeof(double),1,fpo);
  fwrite(&an2,sizeof(double),1,fpo);
  fwrite(&au,sizeof(double),1,fpo);
  fwrite(&aw,sizeof(double),1,fpo);

  fwrite(ag,sizeof(long long),ndr,fpo);

  fclose(fpo);

  /* Deallocate arrays */

  free(ag);
  free(x);
  free(y);
  free(z);

  return;

}

/**********************************************************************/

int main() {

/**********************************************************************/

  int i,j;
  double ran;

  /* Read startup/checkpoint file & initialize variables */

  getcf();

  /* Do ntskip*ntjob passes (particle creations/deletions/moves) */

  for(i=1;i<=ntjob;i++) {
    nt=nt+1;
    for(j=1;j<=n*ntskip;j++) {
      ran=drand48();
      if(ran<pcre) {
	create();               /* Create particle */
      } else if(ran<2.0*pcre) {
	delete();               /* Delete particle */
      } else {
	move();                 /* Move particle */
      }
    }
    means();
  }

  /* Write checkpoint file */

  putcf();

  return 0;

}

/**********************************************************************/
