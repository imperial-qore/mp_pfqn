CC = gcc
LD = gcc
CFLAGS = -g -O4 -fopenmp -I../deps/local/include -I. -I../include/ -I../zpla/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -fopenmp -L../deps/local/lib -lgmp -lm -L. -L../lib/

SRCS = main.c engine.c globals.c image.c setupls.c combs.c btf_decompose.c btf_solve.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o ../gmpla/int_vecsubsum.o
ZPLAOBJS = ../zpla/zpla.o

RM = /bin/rm -f

all: comommod

comommod: ${OBJS} util
	${LD} ${OBJS} ${ZPLAOBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/comommod ${LDFLAGS}

util:  ${UTILOBJS}

include ../Makefile.deps

cleanbase:
	${RM} ${OBJS}

clean:
	${RM} ${OBJS} ${DEPS}
