#**********************************************************************#
#
# Makefile for nptmclj/c
#
# 03-Apr-2010 (MN)
# 04-May-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: anptmclj gnptmclj lnptmclj nptmclj znptmclj

anptmclj: anptmclj.o getval.o
	${CC} ${LFLAGS} -o anptmclj anptmclj.o getval.o

gnptmclj: gnptmclj.o getval.o
	${CC} ${LFLAGS} -o gnptmclj gnptmclj.o getval.o

lnptmclj: lnptmclj.o getval.o
	${CC} ${LFLAGS} -o lnptmclj lnptmclj.o getval.o

nptmclj: nptmclj.c
	${CC} ${CFLAGS} ${LFLAGS} -o nptmclj nptmclj.c

znptmclj: znptmclj.o getval.o
	${CC} ${LFLAGS} -o znptmclj znptmclj.o getval.o

anptmclj.o: anptmclj.c getval.h
	${CC} ${CFLAGS} -c anptmclj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

gnptmclj.o: gnptmclj.c getval.h
	${CC} ${CFLAGS} -c gnptmclj.c

lnptmclj.o: lnptmclj.c getval.h
	${CC} ${CFLAGS} -c lnptmclj.c

znptmclj.o: znptmclj.c getval.h
	${CC} ${CFLAGS} -c znptmclj.c

clean:
	rm -f anptmclj.o getval.o gnptmclj.o lnptmclj.o znptmclj.o

realclean: clean
	rm -f anptmclj gnptmclj lnptmclj nptmclj znptmclj

#**********************************************************************#
