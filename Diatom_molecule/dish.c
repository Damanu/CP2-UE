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
#include <string.h>

/* Global Variables & Parameters */
int n=1,nt=0,ntjob=10000;
double V,l,d12=10;
double roh=0.5,t=100,dt=1,eps=0.000001;
double lamv;

double *vx_1,*vy_1,*vz_1,*vx_2,*vy_2,*vz_2, *vx_1un,*vy_1un,*vz_1un,*vx_2un, *vy_2un, *vz_2un, *vx_1ul,*vy_1ul,*vz_1ul,*vx_2ul, *vy_2ul, *vz_2ul,*x_1,*y_1,*z_1,*x_2,*y_2,*z_2, *x_1ul, *y_1ul, *z_1ul,*x_2ul, *y_2ul, *z_2ul, *x_1un, *y_1un, *z_1un,*x_2un, *y_2un, *z_2un;
double *gx_1,*gy_1,*gz_1,*gx_2,*gy_2,*gz_2;
double *lam,*sigu,*M;

double summe(double* array) {
	int i;
	double sum=0;
	for (i=0;i<n;i++) {
		sum = sum + array[i]; 
	}
	return sum;
}
void init() {
	int i;
	double v_mean;
	/* calculate box measurments and other parameters */
	l=cbrt(n/roh);
	V=l*l*l;
	v_mean=sqrt(3*t);

	for (i=0;i<n;i++) {
		
		//first atom of molecule
		vx_1[i]=drand48()-0.5;
		vy_1[i]=drand48()-0.5;
		vz_1[i]=drand48()-0.5;
	
	/*	vx_1[i]=0;
		vy_1[i]=0;
		vz_1[i]=0;
	*/	
		x_1[i]=drand48()*l;
		y_1[i]=drand48()*l;
		z_1[i]=drand48()*l;
		
/*		x_1[i]=0.0;
		y_1[i]=1.0;
		z_1[i]=0.0;
*/		
//		x_1ul[i]=drand48()*l;
//		y_1ul[i]=drand48()*l;
//		z_1ul[i]=drand48()*l;

		//second atom of molecule
		vx_2[i]=drand48()-0.5;
		vy_2[i]=drand48()-0.5;
		vz_2[i]=drand48()-0.5;
	
	//	vx_2[i]=0;
	//	vy_2[i]=1;
	//	vz_2[i]=0;
	
		x_2[i]=drand48()*l;
		y_2[i]=drand48()*l;
		z_2[i]=drand48()*l;

/*		x_2[i]=0.0;
		y_2[i]=0.0;
		z_2[i]=0.0;
*/
		//recalculate the position of the second atom for right distances between them
//		x_2[i]=d12*(x_2[i]-x_1[i])/(sqrt((x_2[i]-x_1[i])*(x_2[i]-x_1[i])+(y_2[i]-y_1[i])*(y_2[i]-y_1[i])+(z_2[i]-z_1[i])*(z_2[i]-z_1[i])));
	//	y_2[i]=d12*(y_2[i]-y_1[i])/(sqrt((x_2[i]-x_1[i])*(x_2[i]-x_1[i])+(y_2[i]-y_1[i])*(y_2[i]-y_1[i])+(z_2[i]-z_1[i])*(z_2[i]-z_1[i])));
	//	z_2[i]=d12*(z_2[i]-z_1[i])/(sqrt((x_2[i]-x_1[i])*(x_2[i]-x_1[i])+(y_2[i]-y_1[i])*(y_2[i]-y_1[i])+(z_2[i]-z_1[i])*(z_2[i]-z_1[i])));

		//recalculate velocities for the right temperature
/*		vx_1[i]=v_mean*vx_1[i]/summe(vx_1);
		vy_1[i]=v_mean*vy_1[i]/summe(vy_1);
		vz_1[i]=v_mean*vz_1[i]/summe(vz_1);
		
		vx_2[i]=v_mean*vx_2[i]/summe(vx_2);
		vy_2[i]=v_mean*vy_2[i]/summe(vy_2);
		vz_2[i]=v_mean*vz_2[i]/summe(vz_2);*/
//		printf("%10d\t%12.5le\t%12.5le\t%12.5le\t\t%12.5le\t%12.5le\t%12.5le\t\n",i,x_1[i],y_1[i],z_1[i],x_2[i],y_2[i],z_2[i]);
//		printf("%10d\t%12.5le\t%12.5le\t%12.5le\t\t%12.5le\t%12.5le\t%12.5le\t\n",i,vx_1[i],vy_1[i],vz_1[i],vx_2[i],vy_2[i],vz_2[i]);
	}	
}


double absval(double x,double y,double z) {
	int i;
	double val; 
	val=sqrt(x*x+y*y+z*z);
//	printf("%12.5le\n",val);
	return val;
}

void get_lam() {
	int i;
	for(i=0;i<n;i++) {
		sigu[i]=sqrt((x_2ul[i]-x_1ul[i])*(x_2ul[i]-x_1ul[i])+(y_2ul[i]-y_1ul[i])*(y_2ul[i]-y_1ul[i])+(z_2ul[i]-z_1ul[i])*(z_2ul[i]-z_1ul[i]))-d12;
	}
	for(i=0;i<n;i++) {
		M[i]=((x_1ul[i]-x_2ul[i])*(x_1[i]-x_2[i])+(y_1ul[i]-y_2ul[i])*(y_1[i]-y_2[i])+(z_1ul[i]-z_2ul[i])*(z_1[i]-z_2[i])+(x_2ul[i]-x_1ul[i])*(x_2[i]-x_1[i])+(y_2ul[i]-y_1ul[i])*(y_2[i]-y_1[i])+(z_2ul[i]-z_1ul[i])*(z_2[i]-z_1[i]))/(absval(x_2ul[i]-x_1ul[i],y_2ul[i]-y_1ul[i],z_2ul[i]-z_1ul[i])*absval(x_2[i]-x_1[i],y_2[i]-y_1[i],z_2[i]-z_1[i]));
	}
	for(i=0;i<n;i++) {
		lam[i]=sigu[i]/(M[i]*dt*dt);
	}
}
void get_lamv() {
	int i;
	double vx_12,vy_12,vz_12,x_12,y_12,z_12;
	for(i=0;i<n;i++) {
		vx_12=vx_1ul[i]-vx_2ul[i];
		vy_12=vy_1ul[i]-vy_2ul[i];
		vz_12=vz_1ul[i]-vz_2ul[i];
		
		x_12=x_1ul[i]-x_2ul[i];
		y_12=y_1ul[i]-y_2ul[i];
		z_12=z_1ul[i]-z_2ul[i];
		
		lamv=(vx_12*x_12+vy_12*y_12+vz_12*z_12)/(2*(x_12*x_12+y_12*y_12+z_12*z_12));
		
		
		
	}
}
void grad_sig() {
	int i;
	double d12i;
	//calculate gradient of sigma for both atoms
	for(i=0;i<n;i++) {

		//	d12i=sqrt((x_2ul[i]-x_1ul[i])*(x_2ul[i]-x_1ul[i])+(y_2ul[i]-y_1ul[i])*(y_2ul[i]-y_1ul[i])+(z_2ul[i]-z_1ul[i])*(z_2ul[i]-z_1ul[i])); //real distance of atoms
		d12i=absval(x_2ul[i]-x_1ul[i],y_2ul[i]-y_1ul[i],z_2ul[i]-z_1ul[i]);
	//	printf("d12i: %12.5le\n",d12i);	
		gx_1[i]=(x_1ul[i]-x_2ul[i])/d12i;
		gy_1[i]=(y_1ul[i]-y_2ul[i])/d12i;
		gz_1[i]=(z_1ul[i]-z_2ul[i])/d12i;

		gx_2[i]=-gx_1[i];
		gy_2[i]=-gy_1[i];
		gz_2[i]=-gz_1[i];
	}
}

void r_last_new() {
	int i;
	for (i=0;i<n;i++) {
		x_1ul[i]=x_1un[i];	
		y_1ul[i]=y_1un[i];	
		z_1ul[i]=z_1un[i];	
		
		x_2ul[i]=x_2un[i];	
		y_2ul[i]=y_2un[i];	
		z_2ul[i]=z_2un[i];	
	}
}
void v_last_new() {
	int i;
	for (i=0;i<n;i++) {
		vx_1ul[i]=vx_1un[i];	
		vy_1ul[i]=vy_1un[i];	
		vz_1ul[i]=vz_1un[i];	
		
		vx_2ul[i]=vx_2un[i];	
		vy_2ul[i]=vy_2un[i];	
		vz_2ul[i]=vz_2un[i];	
	}
}
void r_last() {
	int i;
	for (i=0;i<n;i++) {
		x_1ul[i]=x_1[i];	
		y_1ul[i]=y_1[i];	
		z_1ul[i]=z_1[i];	
		
		x_2ul[i]=x_2[i];	
		y_2ul[i]=y_2[i];	
		z_2ul[i]=z_2[i];	
	}
}
void v_last() {
	int i;
	for (i=0;i<n;i++) {
		vx_1ul[i]=vx_1[i];	
		vy_1ul[i]=vy_1[i];	
		vz_1ul[i]=vz_1[i];	
		
		vx_2ul[i]=vx_2[i];	
		vy_2ul[i]=vy_2[i];	
		vz_2ul[i]=vz_2[i];	
	}
}
void r_new() {
	int i;
	for (i=0;i<n;i++) {
		x_1[i]=x_1un[i];	
		y_1[i]=y_1un[i];	
		z_1[i]=z_1un[i];	
		
		x_2[i]=x_2un[i];	
		y_2[i]=y_2un[i];	
		z_2[i]=z_2un[i];	
	}
}
void v_new() {
	int i;
	for (i=0;i<n;i++) {
		vx_1[i]=vx_1un[i];	
		vy_1[i]=vy_1un[i];	
		vz_1[i]=vz_1un[i];	
		
		vx_2[i]=vx_2un[i];	
		vy_2[i]=vy_2un[i];	
		vz_2[i]=vz_2un[i];	
	}
}
void move() {
	double d12i,betr;
	int i;
	//velocity verlet integration (force is 0 so there is no change in velocity. True with constraints??)
	for (i=0;i<n;i++) {
		x_1ul[i]=x_1[i]-dt*vx_1[i];
		y_1ul[i]=y_1[i]-dt*vy_1[i];
		z_1ul[i]=z_1[i]-dt*vz_1[i];
		
		x_2ul[i]=x_2[i]-dt*vx_2[i];
		y_2ul[i]=y_2[i]-dt*vy_2[i];
		z_2ul[i]=z_2[i]-dt*vz_2[i];
	}
	v_last();
	do {
		betr=0.0;
		grad_sig();//calculate the gradient of sigma unconstrained new for all molecules
		get_lam();//calculate lambda for all molecules			
		get_lamv();//calculate lambda for all molecules			
		for(i=0;i<n;i++) {
		
			//Rattle: change the coords according to the constraints
			x_1un[i]=x_1ul[i]-dt*dt*lam[i]*gx_1[i];
			y_1un[i]=y_1ul[i]-dt*dt*lam[i]*gy_1[i];
			z_1un[i]=z_1ul[i]-dt*dt*lam[i]*gz_1[i];
			
			x_2un[i]=x_2ul[i]-dt*dt*lam[i]*gx_2[i];
			y_2un[i]=y_2ul[i]-dt*dt*lam[i]*gy_2[i];
			z_2un[i]=z_2ul[i]-dt*dt*lam[i]*gz_2[i];
			
			vx_1un[i]=vx_1ul[i]-dt*dt*lamv*gx_1[i];
			vy_1un[i]=vy_1ul[i]-dt*dt*lamv*gy_1[i];
			vz_1un[i]=vz_1ul[i]-dt*dt*lamv*gz_1[i];
			
			vx_2un[i]=vx_2ul[i]-dt*dt*lamv*gx_2[i];
			vy_2un[i]=vy_2ul[i]-dt*dt*lamv*gy_2[i];
			vz_2un[i]=vz_2ul[i]-dt*dt*lamv*gz_2[i];
				
			//compare the new and the last x, means calculate abs dist between vectors
			betr=betr+sqrt((x_1un[i]-x_1ul[i])*(x_1un[i]-x_1ul[i])+(y_1un[i]-y_1ul[i])*(y_1un[i]-y_1ul[i])+(z_1un[i]-z_1ul[i])*(z_1un[i]-z_1ul[i]));
		}
		r_last_new();//save new x and overwrite the last x
		v_last_new();
		betr=betr/(double)n;
	//	printf("betr: %12.5le\n",betr);	
	//	printf("lam: %12.5le\n",lam[0]);	
	}while(betr>eps);
//	printf("betr: %12.5le\n",betr);	
	r_new();
	v_new();	
}

void means() {
	double d12_mean,u=0;
	int i;
	for(i=0;i<n;i++) {
		u=u+(vx_1[i]*vx_1[i]+vy_1[i]*vy_1[i]+vz_1[i]*vz_1[i]+vx_2[i]*vx_2[i]+vy_2[i]*vy_2[i]+vz_2[i]*vz_2[i])/2.0;
	}
	for(i=0;i<n;i++) {
		d12_mean=d12_mean+absval(x_1[i]-x_2[i],y_1[i]-y_2[i],z_1[i]-z_2[i]);
	}
	d12_mean=d12_mean/(double)n;
	
	u=u/(double)n;
	printf("%10d\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\t%12.5le\n",nt,u,d12_mean,x_1[0],y_1[0],z_1[0],x_2[0],y_2[0],z_2[0],vx_1[0],vy_1[0],vz_1[0],vx_2[0],vy_2[0],vz_2[0]);
}

int main() {
	double abs,r1,r2,r3;
	abs=absval(1.0,1.0,0.0);
	r1=drand48();
	r2=drand48();
	r3=drand48();
//	printf("rand: %12.5le\n",r1);
//	printf("rand: %12.5le\n",r2);
//	printf("rand: %12.5le\n",r3);
	/* allcate arrays */
	vx_1=(double*)malloc(n*sizeof(double));
	vy_1=(double*)malloc(n*sizeof(double));
	vz_1=(double*)malloc(n*sizeof(double));

	vx_2=(double*)malloc(n*sizeof(double));
	vy_2=(double*)malloc(n*sizeof(double));
	vz_2=(double*)malloc(n*sizeof(double));

	vx_1ul=(double*)malloc(n*sizeof(double));
	vy_1ul=(double*)malloc(n*sizeof(double));
	vz_1ul=(double*)malloc(n*sizeof(double));

	vx_2ul=(double*)malloc(n*sizeof(double));
	vy_2ul=(double*)malloc(n*sizeof(double));
	vz_2ul=(double*)malloc(n*sizeof(double));

	vx_1un=(double*)malloc(n*sizeof(double));
	vy_1un=(double*)malloc(n*sizeof(double));
	vz_1un=(double*)malloc(n*sizeof(double));

	vx_2un=(double*)malloc(n*sizeof(double));
	vy_2un=(double*)malloc(n*sizeof(double));
	vz_2un=(double*)malloc(n*sizeof(double));

	x_1=(double*)malloc(n*sizeof(double));
	y_1=(double*)malloc(n*sizeof(double));
	z_1=(double*)malloc(n*sizeof(double));

	x_2=(double*)malloc(n*sizeof(double));
	y_2=(double*)malloc(n*sizeof(double));
	z_2=(double*)malloc(n*sizeof(double));

	x_1ul=(double*)malloc(n*sizeof(double));
	y_1ul=(double*)malloc(n*sizeof(double));
	z_1ul=(double*)malloc(n*sizeof(double));

	x_2ul=(double*)malloc(n*sizeof(double));
	y_2ul=(double*)malloc(n*sizeof(double));
	z_2ul=(double*)malloc(n*sizeof(double));

	x_1un=(double*)malloc(n*sizeof(double));
	y_1un=(double*)malloc(n*sizeof(double));
	z_1un=(double*)malloc(n*sizeof(double));

	x_2un=(double*)malloc(n*sizeof(double));
	y_2un=(double*)malloc(n*sizeof(double));
	z_2un=(double*)malloc(n*sizeof(double));

	gx_1=(double*)malloc(n*sizeof(double));
	gy_1=(double*)malloc(n*sizeof(double));
	gz_1=(double*)malloc(n*sizeof(double));

	gx_2=(double*)malloc(n*sizeof(double));
	gy_2=(double*)malloc(n*sizeof(double));
	gz_2=(double*)malloc(n*sizeof(double));


	lam=(double*)malloc(n*sizeof(double));
	sigu=(double*)malloc(n*sizeof(double));
	M=(double*)malloc(n*sizeof(double));

	init();
	means();
	for(nt=1;nt<=ntjob;nt++) {
		move();
		means();
	}
}



