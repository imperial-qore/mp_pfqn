CC = gcc
LD = gcc
CFLAGS = -g -O2 -I../deps/local/include -I. -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lm

SRCS = main.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../gmpla/int_vec.o ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/readmodel.o

RM = /bin/rm -f
all: automom
automom: ${OBJS} ${UTILOBJS}
	${LD} ${OBJS} ${UTILOBJS} -o ../bin/automom ${LDFLAGS}
include ../Makefile.deps
clean:
	${RM} ${OBJS} ${DEPS}
