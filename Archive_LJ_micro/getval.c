/**********************************************************************
 *
 * File: getval.c
 *
 * Procedures for reading items from stdin/keeping old values
 * on hitting return
 *
 * 03-Apr-2010 (MN)
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getval.h"

#define NL '\n'
#define BSIZE 80
char line[BSIZE];

/**********************************************************************/

void getval_i(int *val) {

/**********************************************************************/

  char *c;

  fgets(line,BSIZE,stdin);
  if((c=(char*)index(line,NL))!=NULL) {
    *c=0;
  }
  if(strlen(line)!=0) {
    *val=atoi(line);
  }

  return;

}

/**********************************************************************/

void getval_d(double *val) {

/**********************************************************************/

  char *c;

  fgets(line,BSIZE,stdin);
  if((c=(char*)index(line,NL))!=NULL) {
    *c=0;
  }
  if(strlen(line)!=0) {
    *val=strtod(line,NULL);
  }

  return;

}

/**********************************************************************/

void getval_s(char *val) {

/**********************************************************************/

  char *c;

  fgets(line,BSIZE,stdin);
  if((c=(char*)index(line,NL))!=NULL) {
    *c=0;
  }
  if(strlen(line)!=0) {
    strcpy(val,line);
  }

  return;

}

/**********************************************************************/
