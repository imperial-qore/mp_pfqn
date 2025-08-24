CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 
LDFLAGS = -L../deps/local/lib -lgmp -lmpfr -lm

SRCS = main.c recal-multi-exact.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/population.o ../util/readmodel.o ../util/printmodel.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o
GMPLAOBJS = ../gmpla/int_vec.o ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o

RM = /bin/rm -f

all: recal

recal: ${OBJS} util
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/recal ${LDFLAGS} 

util: ${UTILOBJS}

# Include automatic dependency generation
include ../Makefile.deps

cleanbase: 
	${RM} ${OBJS}

clean: 
	${RM} ${OBJS} ${DEPS} 