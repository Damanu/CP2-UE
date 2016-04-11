#**********************************************************************#
#
# Makefile for mclj/c
#
# 03-Apr-2010 (MN)
# 19-Mar-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: amclj gmclj lmclj mclj zmclj

amclj: amclj.o getval.o
	${CC} ${LFLAGS} -o amclj amclj.o getval.o

gmclj: gmclj.o getval.o
	${CC} ${LFLAGS} -o gmclj gmclj.o getval.o

lmclj: lmclj.o getval.o
	${CC} ${LFLAGS} -o lmclj lmclj.o getval.o

mclj: mclj.c
	${CC} ${CFLAGS} ${LFLAGS} -o mclj mclj.c

zmclj: zmclj.o getval.o
	${CC} ${LFLAGS} -o zmclj zmclj.o getval.o

amclj.o: amclj.c getval.h
	${CC} ${CFLAGS} -c amclj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

gmclj.o: gmclj.c getval.h
	${CC} ${CFLAGS} -c gmclj.c

lmclj.o: lmclj.c getval.h
	${CC} ${CFLAGS} -c lmclj.c

zmclj.o: zmclj.c getval.h
	${CC} ${CFLAGS} -c zmclj.c

clean:
	rm -f amclj.o getval.o gmclj.o lmclj.o zmclj.o

realclean: clean
	rm -f amclj gmclj lmclj mclj zmclj

#**********************************************************************#
