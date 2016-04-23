#**********************************************************************#
#
# Makefile for gcmclj/c
#
# 08-Apr-2010 (MN)
# 28-Apr-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: agcmclj gcmclj ggcmclj zgcmclj

agcmclj: agcmclj.o getval.o
	${CC} ${LFLAGS} -o agcmclj agcmclj.o getval.o

gcmclj: gcmclj.c
	${CC} ${CFLAGS} ${LFLAGS} -o gcmclj gcmclj.c

ggcmclj: ggcmclj.o getval.o
	${CC} ${LFLAGS} -o ggcmclj ggcmclj.o getval.o

zgcmclj: zgcmclj.o getval.o
	${CC} ${LFLAGS} -o zgcmclj zgcmclj.o getval.o

agmclj.o: agmclj.c getval.h
	${CC} ${CFLAGS} -c agcmclj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

ggmclj.o: ggmclj.c getval.h
	${CC} ${CFLAGS} -c ggcmclj.c

zgmclj.o: zgmclj.c getval.h
	${CC} ${CFLAGS} -c zgcmclj.c

clean:
	rm -f agcmclj.o getval.o ggcmclj.o zgcmclj.o

realclean: clean
	rm -f agcmclj gcmclj ggcmclj zgcmclj

#**********************************************************************#
