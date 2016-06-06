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
double x, v, dt, u, x_0, x_1, v_0, v_1;
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
	
//	x_1=sin(dt);
//	v_1=cos(dt);
}
void move()	{
	double F_0, F_1;

	F_0=-x_0;
	x=x_0+dt*v_0+dt*F_0/2.0;
	F_1=-x;
	v=v_0+dt*(F_0+F_1)/2.0;
	
	v_0=v;
	x_0=x;

	return;	
}
void means()	{
	u=(v_0*v_0+x_0*x_0)/2.0; //total energy
	printf("%10d\t%12.5le\t%12.5le\n",nt,u,x_0);
	return;
}		

int main()	{
	getcf();
	for(nt=0; nt<=ntjob; nt++)	{
		move();
		if(nt%ntprint==0)	{
			means();
		}
	}
}
