CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 
LDFLAGS = -L../deps/local/lib -lgmp -lm -L. -L../lib/ 
#LDFLAGS = -lgmp -lm -lmpfr -lgmpla -L.  -L../lib/
#LDFLAGS = -static -L./gmp-4.1.4/.libs/ -lgmp -lm  -L../lib/

SRCS = blocksolve.c main.c mdecrease.c perfindices.c setupls.c mpq_ludcmp_progress.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_ludcmp.o ../gmpla/mpq_matprint.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_msp.o ../gmpla/mpq_mspget.o ../gmpla/mpq_mspprint.o ../gmpla/mpq_mspset_si.o ../gmpla/mpq_mspset_z.o ../gmpla/mpq_mspvecmul.o ../gmpla/mpq_vec.o ../gmpla/mpq_vecdup.o ../gmpla/mpq_vecprint.o

RM = /bin/rm -f

all: mom

mom: ${OBJS} util
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/mom ${LDFLAGS} 

util:  ${UTILOBJS}

# Include automatic dependency generation
include ../Makefile.deps

cleanbase: 
	${RM} ${OBJS}

clean: 
	${RM} ${OBJS} ${DEPS} 
