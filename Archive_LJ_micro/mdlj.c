/**********************************************************************
 *
 * File: mdlj.c
 *
 * NVE Molecular Dynamics of Lennard-Jonesium
 * Velocity Verlet
 * Pair distribution function g(r)
 * Velocity autocorrelation function C(t)
 *
 * 08-May-2010 (MN)
 * 04-May-2012
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a,b) ((a)<(b)?(a):(b))

/**********************************************************************/

/* Parameters & global variables */

int n,ncor,nt,ntaskip,ntcskip,ntjob,ntorig,ntprint;
int ndr;
long long *ag;
double dr,dt,rho;
double c,c2,c2m1,drm1,r2max,rc,rc2,su0,sw0;
double *fx,*fy,*fz,*vx,*vy,*vz,*x,*y,*z;
double *vxtb,*vytb,*vztb;
double **vxt,**vyt,**vzt;
double ak,ak2,au,aw;
double *acf;

/**********************************************************************/

void getcf() {

/**********************************************************************/

  FILE *fpi;
  int i;

  /* Read startup/checkpoint file */

  fpi=fopen("mdlj_in.dat","r");

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

  /* Allocate arrays */

  ndr=(int)(0.5*cbrt(n/rho)/dr);

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
  vxtb=(double*)malloc((ncor+1)*n*sizeof(double));
  vytb=(double*)malloc((ncor+1)*n*sizeof(double));
  vztb=(double*)malloc((ncor+1)*n*sizeof(double));

  for(i=0;i<=ncor;i++) {
    vxt[i]=vxtb+i*n;
    vyt[i]=vytb+i*n;
    vzt[i]=vztb+i*n;
  }

  /* Positions & velocities */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  fread(vx,sizeof(double),n,fpi);
  fread(vy,sizeof(double),n,fpi);
  fread(vz,sizeof(double),n,fpi);

  /* Accumulated averages/clear accumulators */

  if(nt>0) {
    fread(&ak,sizeof(double),1,fpi);
    fread(&ak2,sizeof(double),1,fpi);
    fread(&au,sizeof(double),1,fpi);
    fread(&aw,sizeof(double),1,fpi);
    fread(ag,sizeof(long long),ndr,fpi);
  } else {
    ak=0.0;
    ak2=0.0;
    au=0.0;
    aw=0.0;
    for(i=0;i<=ndr-1;i++) {
      ag[i]=0;
    }
  }

  /* Ring buffer & velocity correlation function */

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
  su0=2.0*M_PI*rho*n*(4.0/(9.0*pow(rc,9))-4.0/(3.0*pow(rc,3)));
  sw0=2.0*M_PI*rho*n*(24.0/(3.0*pow(rc,3))-48.0/(9.0*pow(rc,9)));

  return;

}

/**********************************************************************/

void force() {

/**********************************************************************/

  int i,j;
  double df,dx,dy,dz,r2,rm6;

  /* Initialize forces */

  for(i=0;i<=n-1;i++) {
    fx[i]=0.0;
    fy[i]=0.0;
    fz[i]=0.0;
  }

  /* All pairs */

  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;        /* Periodic boundary conditions */
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {                   /* In range? */
	rm6=1.0/(r2*r2*r2);
	df=(24.0-48.0*rm6)*rm6/r2;
	fx[i]=fx[i]+df*dx;
	fx[j]=fx[j]-df*dx;           /* Actio=reactio */
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

void halfstep() {

/**********************************************************************/

  int i;

  /* Velocity Verlet */

  /* Calculate new positions & advance velocities by dt/2 */

  for(i=0;i<=n-1;i++) {

    /* New positions (+PBC) */

    x[i]=x[i]+dt*(vx[i]+0.5*dt*fx[i]);
    x[i]=x[i]-(int)(x[i]*c2m1)*c;
    y[i]=y[i]+dt*(vy[i]+0.5*dt*fy[i]);
    y[i]=y[i]-(int)(y[i]*c2m1)*c;
    z[i]=z[i]+dt*(vz[i]+0.5*dt*fz[i]);
    z[i]=z[i]-(int)(z[i]*c2m1)*c;

    /* Intermediate velocities */

    vx[i]=vx[i]+0.5*dt*fx[i];
    vy[i]=vy[i]+0.5*dt*fy[i];
    vz[i]=vz[i]+0.5*dt*fz[i];

  }

  return;

}

/**********************************************************************/

void fullstep() {

/**********************************************************************/

  int i;

  /* Velocity Verlet */

  /* Final velocities */

  for(i=0;i<=n-1;i++) {
    vx[i]=vx[i]+0.5*dt*fx[i];
    vy[i]=vy[i]+0.5*dt*fy[i];
    vz[i]=vz[i]+0.5*dt*fz[i];
  }

  return;

}

/**********************************************************************/

void vacf() {

/**********************************************************************/

  int i,k,nt1,nt2;
  double s;

  /* Store current values in ring buffer (index nt2) */

  nt2=(nt/ntcskip-1)%(ncor+1);

  for(i=0;i<=n-1;i++) {
    vxt[nt2][i]=vx[i];
    vyt[nt2][i]=vy[i];
    vzt[nt2][i]=vz[i];
  }

  /* Correlate with values at previous times (index nt1) */
  /* k = time lag                                        */

  if((nt/ntcskip>ncor)&&((nt/ntcskip-1)%ntorig==0)) {
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
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<r2max) {                   /* In range? */
	if(r2<rc2) {
	  rm6=1.0/(r2*r2*r2);
	  su=su+(4.0*rm6-4.0)*rm6;     /* Potential energy */
	  sw=sw+(24.0-48.0*rm6)*rm6;   /* Virial           */
	}
	k=(int)(sqrt(r2)*drm1);
	if(k<=ndr-1) {
	  ag[k]=ag[k]+1;               /* g(r) */
	}
      }
    }
  }

  /* Accumulate averages */

  ak=ak+sk/n;
  ak2=ak2+(sk/n)*(sk/n);
  au=au+su/n;
  aw=aw+sw/n;

  /* Print control variables: temperature, kinetic, potential & total */
  /* energies, pressure                                               */

  if((ntprint>0)&&((nt/ntaskip)%ntprint==0)) {
    printf(" %10d %12.5le %12.5le %12.5le %12.5le %12.5le\n",nt,\
	   2.0*(sk/n)/3.0,sk/n,su/n,(sk+su)/n,rho*(2.0*sk-sw)/(3*n));
  }

  return;

}

/**********************************************************************/

void putcf() {

/**********************************************************************/

  FILE *fpo;
  int i;

  /* Write checkpoint file */

  fpo=fopen("mdlj_out.dat","w");

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

  fwrite(x,sizeof(double),n,fpo);
  fwrite(y,sizeof(double),n,fpo);
  fwrite(z,sizeof(double),n,fpo);

  fwrite(vx,sizeof(double),n,fpo);
  fwrite(vy,sizeof(double),n,fpo);
  fwrite(vz,sizeof(double),n,fpo);

  fwrite(&ak,sizeof(double),1,fpo);
  fwrite(&ak2,sizeof(double),1,fpo);
  fwrite(&au,sizeof(double),1,fpo);
  fwrite(&aw,sizeof(double),1,fpo);
  fwrite(ag,sizeof(long long),ndr,fpo);

  if(nt/ntcskip>0) {
    for(i=0;i<=min(nt/ntcskip-1,ncor);i++) {
      fwrite(vxt[i],sizeof(double),n,fpo);
      fwrite(vyt[i],sizeof(double),n,fpo);
      fwrite(vzt[i],sizeof(double),n,fpo);
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
  free(x);
  free(y);
  free(z);

  return;

}

/**********************************************************************/

int main() {

/**********************************************************************/

  int i;

  /* Read startup/checkpoint file */

  getcf();

  /* Initialize forces for Velocity Verlet */

  force();

  /* Do ntjob time steps */

  for(i=1;i<=ntjob;i++) {
    nt=nt+1;

    halfstep();           /* Velocity Verlet */
    force();
    fullstep();

    if(nt%ntcskip==0) {   /* Calculate VACF with sampling rate    */
      vacf();             /* 1/(ntcskip*dt), using "time origins" */
    }                     /* ntorig time steps apart              */

    if(nt%ntaskip==0) {   /* Calculate averages every ntaskip */
      means();            /* time steps                       */
    }

  }

  /* Write checkpoint file */

  putcf();

  return 0;

}

/**********************************************************************/
