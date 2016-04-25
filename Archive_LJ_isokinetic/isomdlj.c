/**********************************************************************
 *
 * File: isomdlj.c
 *
 * Isokinetic (NVT) Molecular Dynamics of Lennard-Jonesium
 * Leapfrog
 * Gaussian thermostat
 * Pair distribution function g(r)
 * Velocity autocorrelation function C(t)
 *
 * 20-Nov-2007 (MN)
 * 05-May-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a,b) ((a)<(b)?(a):(b))

/**********************************************************************/

/* Global variables */

int ndr;
int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
long long *ag;
double c,c2,c2m1,drm1,r2max,rc,rc2,su0,sw0,zeta;
double dr,dt,rho,t;
double *fx,*fy,*fz,*vx,*vy,*vz,*wx,*wy,*wz,*x,*y,*z;
double *vxtb,*vytb,*vztb;
double **vxt,**vyt,**vzt;
double ak,au,au2,aw;
double *acf;

/**********************************************************************/

void getcf() {

/**********************************************************************/

  FILE *fpi;
  int i;

  /* Read startup/checkpoint file */

  fpi=fopen("isomdlj_in.dat","r");

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

  /* Allocate arrays */

  ndr=0.5*cbrt(n/rho)/dr;

  acf=(double*)malloc((ncor+1)*sizeof(double));
  ag=(long long*)malloc(ndr*sizeof(long long));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));
  vx=(double*)malloc(n*sizeof(double));
  vy=(double*)malloc(n*sizeof(double));
  vz=(double*)malloc(n*sizeof(double));
  vxt=(double**)malloc((ncor+1)*sizeof(double*));
  vyt=(double**)malloc((ncor+1)*sizeof(double*));
  vzt=(double**)malloc((ncor+1)*sizeof(double*));

  fx=(double*)malloc(n*sizeof(double));
  fy=(double*)malloc(n*sizeof(double));
  fz=(double*)malloc(n*sizeof(double));
  wx=(double*)malloc(n*sizeof(double));
  wy=(double*)malloc(n*sizeof(double));
  wz=(double*)malloc(n*sizeof(double));
  vxtb=(double*)malloc((ncor+1)*n*sizeof(double));
  vytb=(double*)malloc((ncor+1)*n*sizeof(double));
  vztb=(double*)malloc((ncor+1)*n*sizeof(double));

  for(i=0;i<=ncor;i++) {
    vxt[i]=vxtb+i*n;
    vyt[i]=vytb+i*n;
    vzt[i]=vztb+i*n;
  }

  /* Positions & velocities: r(t),v(t-dt/2) */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(vx,sizeof(double),n,fpi);
  fread(vy,sizeof(double),n,fpi);
  fread(vz,sizeof(double),n,fpi);

  /* Read accumulated averages/clear accumulators */

  if(nt>0) {
    fread(&ak,sizeof(double),1,fpi); 
    fread(&au,sizeof(double),1,fpi); 
    fread(&au2,sizeof(double),1,fpi); 
    fread(&aw,sizeof(double),1,fpi);
    fread(ag,sizeof(long long),ndr,fpi);
  } else {
    ak=0.0;
    au=0.0;
    au2=0.0;
    aw=0.0;
    for(i=0;i<ndr-1;i++) {
      ag[i]=0;
    }
  }

  /* Ring buffer & velocity autocorrelation function */

  if(nt/ntcskip>0) {
    for(i=0;i<=min(nt/ntcskip-1,ncor);i++) {
      fread(vxt[i],sizeof(double),n,fpi);
      fread(vyt[i],sizeof(double),n,fpi);
      fread(vzt[i],sizeof(double),n,fpi);
    }
    if(nt/ntcskip>ncor) {
      fread(acf,sizeof(double),ncor+1,fpi);
    }
  } else {
    for(i=0;i<=ncor;i++) {
      acf[i]=0.0;
    }
  }

  fclose(fpi);

  /* Box parameters */

  c=cbrt(n/rho);
  c2=0.5*c;
  c2m1=1.0/c2;
  drm1=1.0/dr;
  r2max=c2*c2;

  /* Cutoff & tail corrections */

  rc=c2;
  rc2=rc*rc;
  su0=2.0*M_PI*rho*n*(4.0/(9.0*pow(rc,9.0))-4.0/(3.0*pow(rc,3.0)));
  sw0=2.0*M_PI*rho*n*(24.0/(3.0*pow(rc,3.0))-48.0/(9.0*pow(rc,9.0)));

  return;

}

/**********************************************************************/

void force() {

/**********************************************************************/

  int i,j;
  double df,dx,dy,dz,r2,rm2,rm6;

  /* Initialize forces */

  for(i=0;i<=n-1;i++) {
    fx[i]=0.0;
    fy[i]=0.0;
    fz[i]=0.0;
  }

  /* Calculate pair forces */

  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(c2m1*dx)*c;       /* Periodic boundary conditions */
      dy=y[j]-y[i];
      dy=dy-(int)(c2m1*dy)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(c2m1*dz)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {                  /* In range? */
	rm2=1.0/r2;
	rm6=rm2*rm2*rm2;
	df=(24.0-48.0*rm6)*rm6*rm2;
 	fx[i]=fx[i]+df*dx;
 	fx[j]=fx[j]-df*dx;          /* Actio=reactio */
 	fy[i]=fy[i]+df*dy;
 	fy[j]=fy[j]-df*dy;
 	fz[i]=fz[i]+df*dz;
 	fz[j]=fz[j]-df*dz;
      }
    }
  }

  return;

}

/**********************************************************************/

void move() {

/**********************************************************************/

  int i;
  double sw2,swf;

  /* Leapfrog + Gaussian thermostat */

  /* Intermediate velocities v(t) */

  swf=0.0;
  sw2=0.0;
  for(i=0;i<=n-1;i++) {
    wx[i]=vx[i]+0.5*dt*fx[i];
    wy[i]=vy[i]+0.5*dt*fy[i];
    wz[i]=vz[i]+0.5*dt*fz[i];
    swf=swf+(wx[i]*fx[i]+wy[i]*fy[i]+wz[i]*fz[i]);
    sw2=sw2+(wx[i]*wx[i]+wy[i]*wy[i]+wz[i]*wz[i]);
  }

  /* Friction coefficient */

  zeta=swf/sw2;

  for(i=0;i<=n-1;i++) {
    vx[i]=vx[i]+dt*(fx[i]-zeta*wx[i]);   /* New velocities v(t+dt/2) */
    vy[i]=vy[i]+dt*(fy[i]-zeta*wy[i]);
    vz[i]=vz[i]+dt*(fz[i]-zeta*wz[i]);
    x[i]=x[i]+dt*vx[i];                  /* New positions & PBC */
    x[i]=x[i]-(int)(c2m1*x[i])*c;
    y[i]=y[i]+dt*vy[i];
    y[i]=y[i]-(int)(c2m1*y[i])*c;
    z[i]=z[i]+dt*vz[i];
    z[i]=z[i]-(int)(c2m1*z[i])*c;
  }

  return;

}

/**********************************************************************/

void vacf() {

/**********************************************************************/

  int i,k,nt1,nt2;
  double s;

  /* Velocity autocorrelation function */

  /* Store current values in ring buffer (index nt2) */

  nt2=(nt/ntcskip-1)%(ncor+1);
  for(i=0;i<=n-1;i++) {
    vxt[nt2][i]=vx[i];
    vyt[nt2][i]=vy[i];
    vzt[nt2][i]=vz[i];
  }

  /* Correlate with values at previous time steps (index nt1)
     k = time lag */

  if((nt/ntcskip>ncor)&&(((nt/ntcskip-1)%ntorig)==0)) {
    for(k=0;k<=ncor;k++) {
      nt1=(nt2-k+ncor+1)%(ncor+1);
      s=0.0;
      for(i=0;i<=n-1;i++) {
	s=s+(vxt[nt1][i]*vxt[nt2][i]+vyt[nt1][i]*vyt[nt2][i]+\
	     vzt[nt1][i]*vzt[nt2][i]);
      }
      acf[k]=acf[k]+s/n;
    }
  }

  return;

}

/**********************************************************************/

void means() {

/**********************************************************************/

  int i,j,k;
  double dx,dy,dz,r2,rm6,sk,su,sw;

  /* Kinetic energy */

  sk=0.0;
  for(i=0;i<=n-1;i++) {
    sk=sk+(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }
  sk=0.5*sk;

  /* Potential energy, virial & g(r) */

  su=su0;
  sw=sw0;

  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {

      dx=x[j]-x[i];
      dx=dx-(int)(c2m1*dx)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(c2m1*dy)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(c2m1*dz)*c;
      r2=dx*dx+dy*dy+dz*dz;

      if(r2<r2max) {                   /* Inscribed sphere */
	k=sqrt(r2)*drm1;
	if(k<ndr) {
	  ag[k]=ag[k]+1;               /* g(r) */
	}
	if(r2<rc2) {                   /* In range? */
	  rm6=1.0/(r2*r2*r2);
	  su=su+(4.0*rm6-4.0)*rm6;     /* Potential energy */
	  sw=sw+(24.0-48.0*rm6)*rm6;   /* Virial           */
	}
      }

    }
  }

  /* Accumulate averages */

  ak=ak+sk/n;
  au=au+su/n;
  au2=au2+(su/n)*(su/n);
  aw=aw+sw/n;

  /* Print control variables: friction coefficient, temperature, */
  /* potential energy & pressure                                 */

  if((ntprint>0)&&((nt/ntaskip)%ntprint==0)) {
    printf("%d %g %g %g %g\n",nt,zeta,sk/(1.5*n),su/n,\
	   rho*(2.0*sk-sw)/(3.0*n));
  }

  return;

}

/**********************************************************************/

void putcf() {

/**********************************************************************/

  FILE *fpo;
  int k;

  /* Write checkpoint file */

  fpo=fopen("isomdlj_out.dat","w");

  fwrite(&n,sizeof(int),1,fpo);
  fwrite(&ncor,sizeof(int),1,fpo);
  fwrite(&nt,sizeof(int),1,fpo);
  fwrite(&ntaskip,sizeof(int),1,fpo);
  fwrite(&ntcskip,sizeof(int),1,fpo);
  fwrite(&ntjob,sizeof(int),1,fpo);
  fwrite(&ntorig,sizeof(int),1,fpo);
  fwrite(&ntprint,sizeof(int),1,fpo);

  fwrite(&dr,sizeof(double),1,fpo);
  fwrite(&dt,sizeof(double),1,fpo);
  fwrite(&rho,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);

  fwrite(x,sizeof(double),n,fpo);
  fwrite(y,sizeof(double),n,fpo);
  fwrite(z,sizeof(double),n,fpo);

  fwrite(vx,sizeof(double),n,fpo);
  fwrite(vy,sizeof(double),n,fpo);
  fwrite(vz,sizeof(double),n,fpo);

  fwrite(&ak,sizeof(double),1,fpo);
  fwrite(&au,sizeof(double),1,fpo);
  fwrite(&au2,sizeof(double),1,fpo);
  fwrite(&aw,sizeof(double),1,fpo);

  fwrite(ag,sizeof(long long),ndr,fpo);

  if(nt/ntcskip>0) {
    for(k=0;k<=min(nt/ntcskip-1,ncor);k++) {
      fwrite(vxt[k],sizeof(double),n,fpo);
      fwrite(vyt[k],sizeof(double),n,fpo);
      fwrite(vzt[k],sizeof(double),n,fpo);
    }
    if(nt/ntcskip>ncor) {
      fwrite(acf,sizeof(double),ncor+1,fpo);
    }
  }

  fclose(fpo);

  /* Deallocate arrays */

  free(acf);
  free(ag);
  free(fx);
  free(fy);
  free(fz);
  free(vx);
  free(vy);
  free(vz);
  free(vxt);
  free(vyt);
  free(vzt);
  free(vxtb);
  free(vytb);
  free(vztb);
  free(wx);
  free(wy);
  free(wz);
  free(x);
  free(y);
  free(z);

  return;

}

/**********************************************************************/

int main() {

/**********************************************************************/

  int it;

  /* Read startup/checkpoint file */

  getcf();

  /* Do ntjob time steps */

  for(it=0;it<=ntjob-1;it++) {
    nt=nt+1;

    force();              /* Leapfrog */
    move();

    if(nt%ntcskip==0) {   /* Calculate VACF with sampling rate    */
      vacf();             /* 1/(ntcskip*dt), using "time origins" */
    }                     /* (ntorig*ntcskip) time steps apart    */

    if(nt%ntaskip==0) {   /* Calculate averages every ntaskip */
      means();            /* time steps                       */
    }

  }

  /* Write checkpoint file */

  putcf();

  return 0;

}

/**********************************************************************/
