/********************************************************************************
*
* File: holf.c
*
* Leapfrog simulation of a 1D harmonic oscillator
* 
* 06-06-2016
*
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getval.c"
/*******************************************************************************/

/*Global variables*/
int ntjob, nt, ntprint;
double x, v, dt, u, x_0, x_1, v_0, v_1, F;
/*******************************************************************************/
void getcf()	{
	FILE *fpi;
	
	/* read startup file*/
	fpi=fopen("ho_in.dat","r");
	
	fread(&dt,sizeof(double),1,fpi);
	fread(&x_0,sizeof(double),1,fpi);
	fread(&v_0,sizeof(double),1,fpi);
	fread(&ntjob,sizeof(int),1,fpi);
	fread(&ntprint,sizeof(int),1,fpi);
	fclose(fpi);
	x_1=sin(dt);
	v_1=cos(dt);
}
void move()	{
	F=-x_1;	//force of harmonic osci
	x=2*x_1-x_0+dt*dt*F;	//verlet algorithm

	v_1=(x-x_0)/(2.0*dt);	//get velocity (for energy)

	x_0=x_1;
	x_1=x;
	return;	
}
void means()	{
	u=(v_1*v_1+x_1*x_1)/2.0; //total energy
	printf("%10d\t%12.5le\t%12.5le\t%12.5le\n",nt,u,x_1,v_1);
	return;
}		

int main()	{
	getcf();
	for (dt=0.001;dt<=1;dt+=0.001) {
		x_1=sin(dt);
		v_1=cos(dt);
		for(nt=0; nt<=1/dt; nt++) {
			move();
		}
		printf("%12.5le\t",dt);
		means();
	}
}
