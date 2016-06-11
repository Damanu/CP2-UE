/**************************************************************************************
*applying the shake algorithm on one diatomic molecule with random initial velocities
*
*11-06-2016
*
*
*
*
*          
*     O 1
*      \
*       \d12
*	 \
*         O 2
*      
*	 
*Emanuel Schwarzhans
***************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Global Variables & Parameters */
int i,n;
double V,l,d12=1;
double roh,t,dt,eps;

/* allcate arrays */
v1x=(double*)malloc(n*sizeof(double));
v1y=(double*)malloc(n*sizeof(double));
v1z=(double*)malloc(n*sizeof(double));

v2x=(double*)malloc(n*sizeof(double));
v2y=(double*)malloc(n*sizeof(double));
v2z=(double*)malloc(n*sizeof(double));

v1xul=(double*)malloc(n*sizeof(double));
v1yul=(double*)malloc(n*sizeof(double));
v1zul=(double*)malloc(n*sizeof(double));

v2xul=(double*)malloc(n*sizeof(double));
v2yul=(double*)malloc(n*sizeof(double));
v2zul=(double*)malloc(n*sizeof(double));

x1=(double*)malloc(n*sizeof(double));
y1=(double*)malloc(n*sizeof(double));
z1=(double*)malloc(n*sizeof(double));

x2=(double*)malloc(n*sizeof(double));
y2=(double*)malloc(n*sizeof(double));
z2=(double*)malloc(n*sizeof(double));

x1ul=(double*)malloc(n*sizeof(double));
y1ul=(double*)malloc(n*sizeof(double));
z1ul=(double*)malloc(n*sizeof(double));

x2ul=(double*)malloc(n*sizeof(double));
y2ul=(double*)malloc(n*sizeof(double));
z2ul=(double*)malloc(n*sizeof(double));

x1un=(double*)malloc(n*sizeof(double));
y1un=(double*)malloc(n*sizeof(double));
z1un=(double*)malloc(n*sizeof(double));

x2un=(double*)malloc(n*sizeof(double));
y2un=(double*)malloc(n*sizeof(double));
z2un=(double*)malloc(n*sizeof(double));

gx1=(double*)malloc(n*sizeof(double));
gy1=(double*)malloc(n*sizeof(double));
gz1=(double*)malloc(n*sizeof(double));

gx2=(double*)malloc(n*sizeof(double));
gy2=(double*)malloc(n*sizeof(double));
gz2=(double*)malloc(n*sizeof(double));

/* calculate box measurments and other parameters */
l=cbrt(n/rho);
V=l*l*l;

v_mean=sqrt(3*t)

void init() {
	for (i=0;i<n;i++) {
		//first atom of molecule
		vx1[i]=drand48()-0.5;
		vy1[i]=drand48()-0.5;
		vz1[i]=drand48()-0.5;
		
		x1[i]=drand48()*l;
		y1[i]=drand48()*l;
		z1[i]=drand48()*l;
		
		x1ul[i]=drand48()*l;
		y1ul[i]=drand48()*l;
		z1ul[i]=drand48()*l;

		//second atom of molecule
		vx2[i]=drand48()-0.5;
		vy2[i]=drand48()-0.5;
		vz2[i]=drand48()-0.5;
	
		x2[i]=drand48()*l;
		y2[i]=drand48()*l;
		z2[i]=drand48()*l;

		//recalculate the position of the second atom for right distances between them
		x2[i]=d12*(x2[i]-x1[i])/(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));
		y2[i]=d12*(y2[i]-y1[i])/(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));
		z2[i]=d12*(z2[i]-z1[i])/(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));

		//recalculate velocities for the right temperature
		vx1[i]=v_mean*vx1[i]/sum(vx1);
		vy1[i]=v_mean*vy1[i]/sum(vy1);
		vz1[i]=v_mean*vz1[i]/sum(vz1);
		
		vx2[i]=v_mean*vx2[i]/sum(vx2);
		vy2[i]=v_mean*vy2[i]/sum(vy2);
		vz2[i]=v_mean*vz2[i]/sum(vz2);
	}	
}
double sum(double* array) {
	double sum=0;
	for (i=0;i<n;i++) {
		sum = sum + array[i]; 
	}
	return sum;
}

void rattle() {
	for(i=0;i<n;i++) {
		sigu[i]=sqrt((x2un[i]-x1un[i])*(x2un[i]-x1un[i])+(y2un[i]-y1un[i])*(y2un[i]-y1un[i])+(z2un[i]-z1un[i])*(z2un[i]-z1un[i]))-d12;
	}
	M=
}

void grad_sig() {
	//calculate gradient of sigma for both atoms
	for(i=0;i<n;i++) {

	d12i=sqrt((x2ul[i]-x1ul[i])*(x2ul[i]-x1ul[i])+(y2ul[i]-y1ul[i])*(y2ul[i]-y1ul[i])+(z2ul[i]-z1ul[i])*(z2ul[i]-z1ul[i])); //real distance of atoms

	gx1[i]=(x1ul[i]-x2ul[i])/d12i;
	gy1[i]=(y1ul[i]-y2ul[i])/d12i;
	gz1[i]=(z1ul[i]-z2ul[i])/d12i;
	
	gx2[i]=-gx1;
	gy2[i]=-gy1;
	gz2[i]=-gz1;
	}
}

void r_last_new() {
	for (i=0;i<n;i++) {
		x1ul[i]=x1un[i];	
		y1ul[i]=y1un[i];	
		z1ul[i]=z1un[i];	
		
		x2ul[i]=x2un[i];	
		y2ul[i]=y2un[i];	
		z2ul[i]=z2un[i];	
	}
}
void r_new() {
	for (i=0;i<n;i++) {
		x1[i]=x1un[i];	
		y1[i]=y1un[i];	
		z1[i]=z1un[i];	
		
		x2[i]=x2un[i];	
		y2[i]=y2un[i];	
		z2[i]=z2un[i];	
	}
}

void move() {
	double d12i,betr;
	
	//velocity verlet integration (force is 0 so there is no change in velocity. True with constraints??)
	for (i=0;i<n;i++) {
		x1ul[i]=2*x1[i]-dt*vx1[i];
		y1ul[i]=2*y1[i]-dt*vy1[i];
		z1ul[i]=2*z1[i]-dt*vz1[i];
		
		x2ul[i]=2*x2[i]-dt*vx2[i];
		y2ul[i]=2*y2[i]-dt*vy2[i];
		z2ul[i]=2*z2[i]-dt*vz2[i];
	}
	do {
		betr=0;
		r_last_new();//save new x and overwrite the last x
		for(i=0;i<n;i++) {
			grad_sig();//calculate the gradient of sigma unconstrained new
			
			//Rattle: change the coords according to the constraints
			x1un[i]=x1ul[i]-dt*dt*lam*gx1[i];
			y1un[i]=y1ul[i]-dt*dt*lam*gy1[i];
			z1un[i]=z1ul[i]-dt*dt*lam*gz1[i];
			
			x2un[i]=x2ul[i]-dt*dt*lam*gx2[i];
			y2un[i]=y2ul[i]-dt*dt*lam*gy2[i];
			z2un[i]=z2ul[i]-dt*dt*lam*gz2[i];
			
			//compare the new and the last x, means calculate abs dist between vectors
			betr=betr+sqrt((x1un[i]-x1ul[i])*(x1un[i]-x1ul[i])+(y1un[i]-y1ul[i])*(y1un[i]-y1ul[i])+(z1un[i]-z1ul[i])*(z1un[i]-z1ul[i]));
		}
		betr=betr/(double)n;
		
	}while(betr>eps);
	r_new();
				
}

int main() {
	

}



