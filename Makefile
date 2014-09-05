#
#
# Make file for compiling HMM code in this directory.
# Author: Tapas Kanungo
# Date: 23 February 1998
# $Id: Makefile,v 1.3 1998/02/23 08:12:35 kanungo Exp kanungo $
# 
#
CFLAGS= -g
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
SRCS= testest.c nrutil.c sample.c hmmutils.c forward.c backward.c baum.c 

final_output=esthmm

all :	$(final_output)
	
esthmm:  nrutil.o hmmutils.o sample.o \
		forward.o backward.o baum.o testest.o nrutil.h hmm.h sample.h 
	 $(CC) -o esthmm  sample.o nrutil.o hmmutils.o \
		forward.o backward.o baum.o testest.o -lm
sample.o: sample.h sample.c
	$(CC) -o sample.o -c sample.c
nrutil.o: nrutil.h nrutil.c
	$(CC) -o nrutil.o -c nrutil.c
hmmutils.o: hmm.h hmmutils.c
	$(CC) -o hmmutils.o -c hmmutils.c
forward.o: hmm.h sample.h nrutil.h forward.c
	$(CC) -o forward.o -c forward.c
backward.o: hmm.h sample.h nrutil.h backward.c
	$(CC) -o backward.o -c backward.c
#baum.o: hmm.h sample.h nrutil.h baum.c
#	$(CC) -o baum.o -c baum.c
clean:
	rm *.o $(final_output)
# DO NOT DELETE THIS LINE -- make depend depends on it.

