CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 
LDFLAGS = -L../deps/local/lib -lgmp -lm

SRCS = main.c recal-multi.c recal-multi-exact.c
TARGET_DIR = ../target

UTILOBJS = ../target/util/population.o ../target/util/readmodel.o ../target/util/printmodel.o ../target/util/nck.o ../target/util/nckinit.o ../target/util/multichoose.o
GMPLAOBJS = ../target/gmpla/int_vec.o ../target/gmpla/int_mat.o ../target/gmpla/int_matmatchrow.o

PROG = recal 
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
	-$(RM) -rf $(TARGET_DIR)
	$(RM) $(OBJS) $(DEPS)
