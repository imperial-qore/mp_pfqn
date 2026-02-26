CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall -Wno-alloc-size-larger-than
LDFLAGS = -L../deps/local/lib -lgmp -lmpfr -lm

SRCS = main.c mvamx.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/population.o ../util/readmodel.o ../util/printmodel.o ../util/freemodel.o ../util/modelclass.o
GMPLAOBJS = ../gmpla/int_vec.o ../gmpla/int_mat.o

PROG = mvamx
RM = /bin/rm -f

all: $(PROG)

$(PROG): $(OBJS) util
	$(LD) $(OBJS) $(UTILOBJS) $(GMPLAOBJS) $(LDFLAGS) -o ../bin/$(PROG)

util: $(UTILOBJS)

# Include automatic dependency generation
include ../Makefile.deps

cleanbase:
	$(RM) $(OBJS)

clean:
	$(RM) $(OBJS) $(DEPS)
