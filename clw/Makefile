CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lmpfr -lm

SRCS = main.c pfqn_clw.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/population.o ../util/readmodel.o ../util/printmodel.o ../util/freemodel.o ../util/modelclass.o
GMPLAOBJS = ../gmpla/int_vec.o ../gmpla/int_mat.o

RM = /bin/rm -f

all: clw

clw: $(OBJS) util
	$(LD) $(OBJS) $(UTILOBJS) $(GMPLAOBJS) $(LDFLAGS) -o ../bin/clw

util: $(UTILOBJS)

# Include automatic dependency generation
include ../Makefile.deps

clean:
	$(RM) $(OBJS) $(DEPS)
