#**********************************************************************#
#
# Makefile for isomdlj/c
#
# 28-May-2010
# 05-May-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: aisomdlj lisomdlj isomdlj zisomdlj

aisomdlj: aisomdlj.o getval.o
	${CC} ${LFLAGS} -o aisomdlj aisomdlj.o getval.o

lisomdlj: lisomdlj.o getval.o
	${CC} ${LFLAGS} -o lisomdlj lisomdlj.o getval.o

isomdlj: isomdlj.c
	${CC} ${CFLAGS} ${LFLAGS} -o isomdlj isomdlj.c

zisomdlj: zisomdlj.o getval.o
	${CC} ${LFLAGS} -o zisomdlj zisomdlj.o getval.o

aisomdlj.o: aisomdlj.c getval.h
	${CC} ${CFLAGS} -c aisomdlj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

lisomdlj.o: lisomdlj.c getval.h
	${CC} ${CFLAGS} -c lisomdlj.c

zisomdlj.o: zisomdlj.c getval.h
	${CC} ${CFLAGS} -c zisomdlj.c

clean:
	rm -f aisomdlj.o getval.o lisomdlj.o zisomdlj.o

realclean: clean
	rm -f aisomdlj lisomdlj isomdlj zisomdlj

#**********************************************************************#
