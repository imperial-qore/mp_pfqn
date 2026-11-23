CC = gcc
LD = gcc
CFLAGS = -g -O2 -I../deps/local/include -I. -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lm -L. -L../lib/

SRCS = main.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../gmpla/int_vec.o ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/mpq_ludcmp.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_vec.o

RM = /bin/rm -f
all: gcomom
gcomom: ${OBJS} ${UTILOBJS}
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/gcomom ${LDFLAGS}
include ../Makefile.deps
clean:
	${RM} ${OBJS} ${DEPS}
