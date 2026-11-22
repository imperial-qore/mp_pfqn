CC = gcc
LD = gcc
CFLAGS = -g -O4 -fopenmp -I../deps/local/include -I. -I../include/ -I../zpla/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -fopenmp -L../deps/local/lib -lgmp -lm -L. -L../lib/

SRCS = blocksolve.c main.c engine.c globals.c mdecrease.c exportres.c setupls.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o
ZPLAOBJS = ../zpla/zpla.o

RM = /bin/rm -f

all: mommod

mommod: ${OBJS} util
# -fopenmp on the link line itself: the top-level Makefile passes LDFLAGS on
# the command line, which replaces the assignment above and with it -fopenmp
	${LD} ${OBJS} ${ZPLAOBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/mommod ${LDFLAGS} -fopenmp

util:  ${UTILOBJS}

include ../Makefile.deps

cleanbase:
	${RM} ${OBJS}

clean:
	${RM} ${OBJS} ${DEPS}
