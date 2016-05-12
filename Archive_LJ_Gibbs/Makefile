#**********************************************************************#
#
# Makefile for gemclj/c
#
# 17-Apr-2010 (MN)
# 28-Apr-2012
#
#**********************************************************************#

CC=gcc
CFLAGS=-O
LFLAGS=-lm

all: agemclj gemclj ggemclj zgemclj

agemclj: agemclj.o getval.o
	${CC} ${LFLAGS} -o agemclj agemclj.o getval.o

gemclj: gemclj.c
	${CC} ${CFLAGS} ${LFLAGS} -o gemclj gemclj.c

ggemclj: ggemclj.o getval.o
	${CC} ${LFLAGS} -o ggemclj ggemclj.o getval.o

zgemclj: zgemclj.o getval.o
	${CC} ${LFLAGS} -o zgemclj zgemclj.o getval.o

agemclj.o: agemclj.c getval.h
	${CC} ${CFLAGS} -c agemclj.c

getval.o: getval.c getval.h
	${CC} ${CFLAGS} -c getval.c

ggemclj.o: ggemclj.c getval.h
	${CC} ${CFLAGS} -c ggemclj.c

zgemclj.o: zgemclj.c getval.h
	${CC} ${CFLAGS} -c zgemclj.c

clean:
	rm -f agemclj.o getval.o ggemclj.o zgemclj.o

realclean: clean
	rm -f agemclj gemclj ggemclj zgemclj

#**********************************************************************#
