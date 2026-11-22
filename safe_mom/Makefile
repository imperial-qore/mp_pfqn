CC = gcc
LD = gcc
CFLAGS = -g -O2 -I../deps/local/include -I. -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lm -L. -L../lib/
SRCS = safe_mom.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)
UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o ../util/safe_tiers.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_vec.o
RM = /bin/rm -f
all: safe_mom
safe_mom: ${OBJS} ${UTILOBJS}
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/safe_mom ${LDFLAGS}
include ../Makefile.deps
clean:
	${RM} ${OBJS} ${DEPS}
