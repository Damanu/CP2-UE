/*************************************************************
*
*Initialisation for holf.c
*
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "getval.c"

#define BSIZE 80

int main()	{
	FILE *fpo;
	char fname[BSIZE];
	int ntjob, ntprint;
	double x_0, v_0, dt;
	
	/*user input*/
	printf("\tdt=");
	getval_d(&dt);
	printf("\tx_0=");
	getval_d(&x_0);
	printf("\tv_0=");
	getval_d(&v_0);
	printf("\tntjob=");
	getval_i(&ntjob);
	printf("\tntprint=");
	getval_i(&ntprint);

	strcpy(fname,"ho_in.dat");
	
/* Write startup file */
	fpo=fopen(fname,"w");
	fwrite(&dt,sizeof(double),1,fpo);
	fwrite(&x_0,sizeof(double),1,fpo);
	fwrite(&v_0,sizeof(double),1,fpo);
	fwrite(&ntjob,sizeof(int),1,fpo);
	fwrite(&ntprint,sizeof(int),1,fpo);
	
	fclose(fpo);	
	return 0;
	
}
