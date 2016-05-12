/**********************************************************************
 *
 * File: gemclj.c
 *
 * Gibbs Ensemble (NVT) Monte Carlo of Lennard-Jonesium
 *
 * 01-May-2010 (MN)
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

const double r2min=0.25*0.25;

int n,n1,nt,ntjob,ntprint,ntskip;
int n2,naccp1,naccp2,nacct,naccv,ntryp1,ntryp2,ntryv,ntryt1,ntryt2,\
  nacct,naccv;
double disp1,disp2,dv,ptr,pvm,t,v,v1;
double alnmax,c1,c2,rc1,rc2,smu1,smu2,su01,su02,smu1,smu2,v2;
double *x,*y,*z;
double *uattb,*uattn,*urepb,*urepn;
double **uatt,**urep;
long long accrp1,accrp2,accrt,accrv,atryp1,atryp2,atryt1,atryt2,atryv,\
  an1,an2;
double amu1,ap1,arho1,au1,av1,amu2,ap2,arho2,au2,av2;

/**********************************************************************/

void getcf() {

/**********************************************************************/

  FILE *fpi;
  unsigned short seed[3];
  int i;

  /* Maximum argument for exponential */

  alnmax=log(0.1*DBL_MAX);

  /* Read startup/checkpoint file */

  fpi=fopen("gemclj_in.dat","r");

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

  uatt=(double**)malloc(n*sizeof(double*));
  urep=(double**)malloc(n*sizeof(double*));
  uattb=(double*)malloc(n*n*sizeof(double));
  urepb=(double*)malloc(n*n*sizeof(double));
  uattn=(double*)malloc(n*sizeof(double));
  urepn=(double*)malloc(n*sizeof(double));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  for(i=0;i<=n-1;i++) {
    urep[i]=urepb+i*n;
    uatt[i]=uattb+i*n;
  }

  /* Positions */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  /* Read accumulated averages/clear accumulators */

  if(nt>0) {
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
  } else {
    accrp1=0;
    accrp2=0;
    accrt=0;
    accrv=0;
    atryp1=0;
    atryp2=0;
    atryt1=0;
    atryt2=0;
    atryv=0;
    amu1=0.0;
    an1=0;
    ap1=0.0;
    arho1=0.0;
    au1=0.0;
    av1=0.0;
    amu2=0.0;
    an2=0;
    ap2=0.0;
    arho2=0.0;
    au2=0.0;
    av2=0.0;
  }

  fclose(fpi);

  /* Box parameters */

  v2=v-v1;
  c1=cbrt(v1);
  c2=cbrt(v2);

  /* Cutoff & tail corrections */

  n2=n-n1;
  rc1=0.5*c1;
  rc2=0.5*c2;
  su01=2.0*M_PI*n1*n1/v1*(4.0/(9.0*pow(rc1,9))-4.0/(3.0*pow(rc1,3)));
  su02=2.0*M_PI*n2*n2/v2*(4.0/(9.0*pow(rc2,9))-4.0/(3.0*pow(rc2,3)));

  /* Initialize RNG & acceptance counters */

  seed48(seed);

  ntryp1=0;
  ntryp2=0;
  ntryv=0;
  ntryt1=0;
  ntryt2=0;
  naccp1=0;
  naccp2=0;
  naccv=0;
  nacct=0;

  smu1=0.0;
  smu2=0.0;

  return;

}

/**********************************************************************/

void uinit() {

/**********************************************************************/

  int i,ibox,imin,imax,j;
  double c,dx,dy,dz,f2c,r2,rcrc,rm6;

  /* Initialize pair interaction matrix
   *
   * urep = 4/r^12, uatt = -4/r^6
   *
   */

  for(i=0;i<=n-1;i++) {
    for(j=0;j<=n-1;j++) {
      urep[i][j]=0.0;
      uatt[i][j]=0.0;
    }
  }

  /* Do box 1 & 2 */

  for(ibox=1;ibox<=2;ibox++) {

    if(ibox==1) {
      imin=0;
      imax=n1-1;
      c=c1;
      rcrc=rc1*rc1;
    } else {
      imin=n1;
      imax=n-1;
      c=c2;
      rcrc=rc2*rc2;
    }
    f2c=2.0/c;

    for(i=imin;i<=imax-1;i++) {
      for(j=i+1;j<=imax;j++) {
	dx=x[j]-x[i];
	dx=dx-(int)(f2c*dx)*c;
	dy=y[j]-y[i];
	dy=dy-(int)(f2c*dy)*c;
	dz=z[j]-z[i];
	dz=dz-(int)(f2c*dz)*c;
	r2=dx*dx+dy*dy+dz*dz;
	if(r2<rcrc) {
	  rm6=1.0/(r2*r2*r2);
	  urep[i][j]=4.0*rm6*rm6;
	  uatt[i][j]=-4.0*rm6;
	  urep[j][i]=urep[i][j];
	  uatt[j][i]=uatt[i][j];
	}
      }
    }

  }

  return;

}

/**********************************************************************/

void move() {

/**********************************************************************/

  int accept;
  int i,j,jmax,jmin;
  double c,disp,dx,dy,dz,f2c,r2,rcrc,rm6,su,xin,yin,zin;

  /* Random particle */

  i=min((int)(drand48()*n),n-1);

  /* Which box? */

  if(i<=n1-1) {
    ntryp1=ntryp1+1;   /* Box 1 */
    jmin=0;
    jmax=n1-1;
    c=c1;
    disp=disp1;
    rcrc=rc1*rc1;
  } else {
    ntryp2=ntryp2+1;   /* Box 2 */
    jmin=n1;
    jmax=n-1;
    c=c2;
    disp=disp2;
    rcrc=rc2*rc2;
  }
  f2c=2.0/c;

  /* Trial move (displacement) */

  xin=x[i]+disp*(drand48()-0.5);
  xin=xin-(int)(f2c*xin)*c;
  yin=y[i]+disp*(drand48()-0.5);
  yin=yin-(int)(f2c*yin)*c;
  zin=z[i]+disp*(drand48()-0.5);
  zin=zin-(int)(f2c*zin)*c;

  /* Energy difference */

  su=0.0;
  for(j=jmin;j<=jmax;j++) {
    if(j==i) {
      urepn[j]=0.0;
      uattn[j]=0.0;
    } else {
      dx=x[j]-xin;
      dx=dx-(int)(f2c*dx)*c;
      dy=y[j]-yin;
      dy=dy-(int)(f2c*dy)*c;
      dz=z[j]-zin;
      dz=dz-(int)(f2c*dz)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rcrc) {
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

  /* Accept move & update interaction matrix */

  if(accept) {

    if(i<=n1-1) {
      naccp1=naccp1+1;
    } else {
      naccp2=naccp2+1;
    }
    x[i]=xin;
    y[i]=yin;
    z[i]=zin;

    for(j=jmin;j<=jmax;j++) {
      urep[i][j]=urepn[j];
      urep[j][i]=urepn[j];
      uatt[i][j]=uattn[j];
      uatt[j][i]=uattn[j];
    }

  }

  return;

}

/**********************************************************************/

void ptrans() {

/**********************************************************************/

  int accept;
  int i,iold,j,jmax,jmin;
  double c,dx,dy,dz,f2c,r2,rcrc,rm6,su,xin,yin,zin;

  /* Direction of transfer */

  if(drand48()<0.5) {
    ntryt1=ntryt1+1;                   /* Box 2 -> box 1 (n1 -> n1+1) */
    if(n1>=n) {
      i=n;                             /* Box 1 full */
    } else {
      i=min(n1+(int)(drand48()*n2),n-1);   /* Delete i from box 2 */
      iold=n1;                             /* Add as n1 to box 1 */
    }
  } else {
    ntryt2=ntryt2+1;                   /* Box 1 -> box 2 (n1 -> n1-1) */
    if(n1<=0) {
      i=-1;                            /* Box 1 full */
    } else {
      i=min((int)(drand48()*n1),n1-1);   /* Delete i from box 1 */
      iold=n1-1;                         /* Add as n1-1 to box 2 */
    }
  }

  /* Test particle energy */

  if(i>=n1) {
    jmin=0;      /* Box 1 */
    jmax=n1-1;
    c=c1;
    rcrc=rc1*rc1;
    su=2.0*M_PI*(2*n1+1)/v1*(4.0/(9.0*pow(rc1,9))-4.0/(3.0*pow(rc1,3)));
  } else {
    jmin=n1;     /* Box 2 */
    jmax=n-1;
    c=c2;
    rcrc=rc2*rc2;
    su=2.0*M_PI*(2*n2+1)/v2*(4.0/(9.0*pow(rc2,9))-4.0/(3.0*pow(rc2,3)));
  }
  f2c=2.0/c;

  xin=c*(drand48()-0.5);   /* Trial position */
  yin=c*(drand48()-0.5);
  zin=c*(drand48()-0.5);

  for(j=jmin;j<=jmax;j++) {   /* Test particle energy */
    dx=x[j]-xin;
    dx=dx-(int)(f2c*dx)*c;
    dy=y[j]-yin;
    dy=dy-(int)(f2c*dy)*c;
    dz=z[j]-zin;
    dz=dz-(int)(f2c*dz)*c;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<rcrc) {
      if(r2<r2min) {     /* Overlap? */
	return;
      } else {
	rm6=1.0/(r2*r2*r2);
	urepn[j]=4.0*rm6*rm6;
	uattn[j]=-4.0*rm6;
	su=su+(urepn[j]+uattn[j]);
      }
    } else {
      urepn[j]=0.0;
      uattn[j]=0.0;
    }
  }

  /* Chemical potential */

  if(su/t<alnmax) {
    if(i>=n1) {                               /* Box 1 */
      if(n1<n) {
	smu1=smu1+exp(-su/t)*v1/(n1+1);
      } else {
	smu1=smu1+2.0*exp(-su/t)*v1/(n1+1);
      }
    } else {                                  /* Box 2 */
      if(n2<n) {
	smu2=smu2+exp(-su/t)*v2/(n2+1);
      } else {
	smu2=smu2+2.0*exp(-su/t)*v2/(n2+1);
      }
    }
  }

  /* Virtual transfer only? */

  if((i<0)||(i>=n)) {
    return;
  }

  /* Complete trial "energy" */

  if(i>=n1) {     /* Delete i from box 2 */
    su=su+\
      2.0*M_PI*(1-2*n2)/v2*(4.0/(9.0*pow(rc2,9))-4.0/(3.0*pow(rc2,3)))+\
      t*log((n1+1)*v2/(n2*v1));
    for(j=n1;j<=n-1;j++) {
      su=su-(urep[i][j]+uatt[i][j]);
    } 
  } else {        /* Delete i from box 1 */
    su=su+\
      2.0*M_PI*(1-2*n1)/v1*(4.0/(9.0*pow(rc1,9))-4.0/(3.0*pow(rc1,3)))+\
      t*log((n2+1)*v1/(n1*v2));
    for(j=0;j<=n1-1;j++) {
      su=su-(urep[i][j]+uatt[i][j]);
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

  /* Accept transfer */

  if(accept) {
    nacct=nacct+1;

    x[i]=x[iold];                   /* Swap positions */
    y[i]=y[iold];
    z[i]=z[iold];
    x[iold]=xin;
    y[iold]=yin;
    z[iold]=zin;

    urepn[i]=urepn[iold];          /* Update pair interaction matrix */
    uattn[i]=uattn[iold];
    for(j=0;j<=n-1;j++) {
      urep[i][j]=urep[iold][j];
      urep[j][i]=urep[i][j];
      urep[iold][j]=urepn[j];
      urep[j][iold]=urep[iold][j];
      uatt[i][j]=uatt[iold][j];
      uatt[j][i]=uatt[i][j];
      uatt[iold][j]=uattn[j];
      uatt[j][iold]=uatt[iold][j];
    }
    urep[i][i]=0.0;
    urep[iold][iold]=0.0;
    uatt[i][i]=0.0;
    uatt[iold][iold]=0.0;

    if(i>=n1) {                /* Particle numbers & tail corrections */
      n1=n1+1;
    } else {
      n1=n1-1;
    }
    n2=n-n1;
    su01=2.0*M_PI*n1*n1/v1*(4.0/(9.0*pow(rc1,9))-4.0/(3.0*pow(rc1,3)));
    su02=2.0*M_PI*n2*n2/v2*(4.0/(9.0*pow(rc2,9))-4.0/(3.0*pow(rc2,3)));

  }
  
  return;

}

/**********************************************************************/

void vmove() {

/**********************************************************************/

  int accept;
  int i,j;
  double fact,fatt,frep,rcn,su,su0n,v1n,v2n;

  ntryv=ntryv+1;

  /* Volume change */

  v1n=v1+(drand48()-0.5)*dv;
  if((v1n<=0.0)||(v1n>=v)) {
    return;
  }
  v2n=v-v1n;

  /* "Energy" difference: Box 1 */

  rcn=0.5*cbrt(v1n);   /* ! New cutoff & tail correction */
  su0n=2.0*M_PI*n1*n1/v1n*(4.0/(9.0*pow(rcn,9))-4.0/(3.0*pow(rcn,3)));

  frep=pow(v1/v1n,4)-1.0;
  fatt=pow(v1/v1n,2)-1.0;
  su=su0n-su01;
  for(i=0;i<=n1-2;i++) {
    for(j=i+1;j<=n1-1;j++) {
      su=su+(frep*urep[i][j]+fatt*uatt[i][j]);
    }
  }

  su=su-t*n1*log(v1n/v1);

  /* Box 2 */

  rcn=0.5*cbrt(v2n);   /* ! New cutoff & tail correction */
  su0n=2.0*M_PI*n2*n2/v2n*(4.0/(9.0*pow(rcn,9))-4.0/(3.0*pow(rcn,3)));

  frep=pow(v2/v2n,4)-1.0;
  fatt=pow(v2/v2n,2)-1.0;
  su=su+(su0n-su02);
  for(i=n1;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      su=su+(frep*urep[i][j]+fatt*uatt[i][j]);
    }
  }

  su=su-t*n2*log(v2n/v2);

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

    fact=cbrt(v1n/v1);      /* Box 1 */
    for(i=0;i<=n1-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }

    fact=cbrt(v2n/v2);      /* Box 2 */
    for(i=n1;i<=n-1;i++) {
      x[i]=fact*x[i];
      y[i]=fact*y[i];
      z[i]=fact*z[i];
    }

    /* Parameters */

    v1=v1n;       /* Box 1 */
    c1=cbrt(v1);
    v2=v2n;       /* Box 2 */
    c2=cbrt(v2);

    /* Cutoff & tail corrections */

    rc1=0.5*c1;     /* Box 1 */
    su01=2.0*M_PI*n1*n1/v1*(4.0/(9.0*pow(rc1,9))-4.0/(3.0*pow(rc1,3)));
    rc2=0.5*c2;     /* Box 2 */
    su02=2.0*M_PI*n2*n2/v2*(4.0/(9.0*pow(rc2,9))-4.0/(3.0*pow(rc2,3)));

    /* Pair interaction matrix */

    uinit();

  }

  return;

}

/**********************************************************************/

void means() {

/**********************************************************************/

  int i,j;
  double c,dx,dy,dz,f2c,r2,rcrc,rho1,rho2,rm6,su1,su2,sw1,sw2;

  /* Potential energy & virial */

  su1=su01;
  su2=su02;
  sw1=2.0*M_PI*n1*n1/v1*(24.0/(3.0*pow(rc1,3))-48.0/(9.0*pow(rc1,9)));
  sw2=2.0*M_PI*n2*n2/v2*(24.0/(3.0*pow(rc2,3))-48.0/(9.0*pow(rc2,9)));

  /* Box 1 */

  c=c1;
  f2c=2.0/c;
  rcrc=rc1*rc1;
  for(i=0;i<=n1-2;i++) {
    for(j=i+1;j<=n1-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(f2c*dx)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(f2c*dy)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(f2c*dz)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rcrc) {
	su1=su1+(urep[i][j]+uatt[i][j]);
	rm6=1.0/(r2*r2*r2);
	sw1=sw1+(24.0-48.0*rm6)*rm6;
      }
    }
  }
  su1=su1/max(n1,1);

  /* Box 2 */

  c=c2;
  f2c=2.0/c;
  rcrc=rc2*rc2;
  for(i=n1;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(f2c*dx)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(f2c*dy)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(f2c*dz)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rcrc) {
	su2=su2+(urep[i][j]+uatt[i][j]);
	rm6=1.0/(r2*r2*r2);
	sw2=sw2+(24.0-48.0*rm6)*rm6;
      }
    }
  }
  su2=su2/max(n2,1);

  /* Print control variables */

  rho1=n1/v1;
  rho2=n2/v2;

  if((ntprint>0)&&(nt%ntprint==0)) {
    printf("%d\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t",nt,\
	   (double)naccp1/max(ntryp1,1),\
	   (double)naccp2/max(ntryp2,1),\
	   (double)nacct/max(ntryt1+ntryt2,1),\
	   (double)naccv/max(ntryv,1));
    printf("%12.5le\t%12.5le\t%12.5le\t%12.5le\n",\
	   su1,su2,rho1,rho2);
  }

  /* Accumulate averages */

  accrp1=accrp1+naccp1;
  accrp2=accrp2+naccp2;
  accrv=accrv+naccv;
  accrt=accrt+nacct;
  atryp1=atryp1+ntryp1;
  atryp2=atryp2+ntryp2;
  atryt1=atryt1+ntryt1;
  atryt2=atryt2+ntryt2;
  atryv=atryv+ntryv;
  an1=an1+n1;
  an2=an2+n2;

  amu1=amu1+smu1;
  ap1=ap1+(t*rho1-sw1/(3.0*v1));
  arho1=arho1+rho1;
  au1=au1+su1;
  av1=av1+v1;

  amu2=amu2+smu2;
  ap2=ap2+(t*rho2-sw2/(3.0*v2));
  arho2=arho2+rho2;
  au2=au2+su2;
  av2=av2+v2;

  /* Clear acceptance counters & accumulators */

  ntryp1=0;
  ntryp2=0;
  ntryt1=0;
  ntryt2=0;
  ntryv=0;
  naccp1=0;
  naccp2=0;
  nacct=0;
  naccv=0;

  smu1=0.0;
  smu2=0.0;

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

  fpo=fopen("gemclj_out.dat","w");

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

  fwrite(&accrp1,sizeof(long long),1,fpo);
  fwrite(&accrp2,sizeof(long long),1,fpo);
  fwrite(&accrt,sizeof(long long),1,fpo);
  fwrite(&accrv,sizeof(long long),1,fpo);
  fwrite(&atryp1,sizeof(long long),1,fpo);
  fwrite(&atryp2,sizeof(long long),1,fpo);
  fwrite(&atryt1,sizeof(long long),1,fpo);
  fwrite(&atryt2,sizeof(long long),1,fpo);
  fwrite(&atryv,sizeof(long long),1,fpo);
  fwrite(&amu1,sizeof(double),1,fpo);
  fwrite(&an1,sizeof(long long),1,fpo);
  fwrite(&ap1,sizeof(double),1,fpo);
  fwrite(&arho1,sizeof(double),1,fpo);
  fwrite(&au1,sizeof(double),1,fpo);
  fwrite(&av1,sizeof(double),1,fpo);
  fwrite(&amu2,sizeof(double),1,fpo);
  fwrite(&an2,sizeof(long long),1,fpo);
  fwrite(&ap2,sizeof(double),1,fpo);
  fwrite(&arho2,sizeof(double),1,fpo);
  fwrite(&au2,sizeof(double),1,fpo);
  fwrite(&av2,sizeof(double),1,fpo);

  fclose(fpo);

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
  double ran;

  /* Read startup/checkpoint file & initialize */

  getcf();
  uinit();

  /* Do ntskip*ntjob passes (transfers/volume changes/displacements) */

 // printf("nt\tacc1\tacc2\tacct\taccv\tsu1\tsu2\trho1\trho2\n");
  for(i=1;i<=ntjob;i++) {
    nt=nt+1;
    for(j=1;j<=ntskip*n;j++) {
      ran=drand48();
      if(ran<ptr) {
	ptrans();                   /* Partice transfer */
      } else if(ran<(ptr+pvm)) {
	vmove();                    /* Volume change */
      } else {
	move();                     /* Particle displacement */
      }
    }
    means();                        /* Analyze configuration */
  }

  /* Write checkpoint file */

  putcf();

  return 0;

}

/**********************************************************************/
