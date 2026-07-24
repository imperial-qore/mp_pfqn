CC = gcc
LD = gcc
CFLAGS = -g -O2 -I../deps/local/include -I. -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lm -L. -L../lib/

SRCS = base.c setup1.c main.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_vec.o ../gmpla/mpq_ludcmp.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_vec.o ../gmpla/mpq_matvecmul.o ../gmpla/mpq_mattransmul.o ../gmpla/mpq_matprint.o

RM = /bin/rm -f
all: gmom
gmom: ${OBJS} ${UTILOBJS}
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/gmom ${LDFLAGS}
include ../Makefile.deps
clean:
	${RM} ${OBJS} ${DEPS}
