CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../fpla/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lmpfr -lgmp -lm -L. -L../lib/

SRCS = blocksolve.c main.c mdecrease.c perfindices.c setupls.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o
FPLAOBJS = ../fpla/fpla.o

RM = /bin/rm -f

all: momf

momf: ${OBJS} util
	${LD} ${OBJS} ${FPLAOBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/momf ${LDFLAGS}

util:  ${UTILOBJS}

include ../Makefile.deps

cleanbase:
	${RM} ${OBJS}

clean:
	${RM} ${OBJS} ${DEPS}
