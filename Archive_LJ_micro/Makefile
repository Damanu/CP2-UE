#**********************************************************************#
#
# Makefile for mdlj/c
#
# 07-May-2010 (MN)
# 04-May-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: amdlj gmdlj lmdlj mdlj zmdlj

amdlj: amdlj.o getval.o
	${CC} ${LFLAGS} -o amdlj amdlj.o getval.o

gmdlj: gmdlj.o getval.o
	${CC} ${LFLAGS} -o gmdlj gmdlj.o getval.o

lmdlj: lmdlj.o getval.o
	${CC} ${LFLAGS} -o lmdlj lmdlj.o getval.o

mdlj: mdlj.c
	${CC} ${CFLAGS} ${LFLAGS} -o mdlj mdlj.c

zmdlj: zmdlj.o getval.o
	${CC} ${LFLAGS} -o zmdlj zmdlj.o getval.o

amdlj.o: amdlj.c getval.h
	${CC} ${CFLAGS} -c amdlj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

gmdlj.o: gmdlj.c getval.h
	${CC} ${CFLAGS} -c gmdlj.c

lmdlj.o: lmdlj.c getval.h
	${CC} ${CFLAGS} -c lmdlj.c

zmdlj.o: zmdlj.c getval.h
	${CC} ${CFLAGS} -c zmdlj.c

clean:
	rm -f amdlj.o getval.o gmdlj.o lmdlj.o zmdlj.o

realclean: clean
	rm -f amdlj gmdlj lmdlj mdlj zmdlj

#**********************************************************************#
