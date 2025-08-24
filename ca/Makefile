CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../mom/ -I../gmpla/ -I../util/ -Wall -Wno-alloc-size-larger-than 
LDFLAGS = -L../deps/local/lib -lgmp -lmpfr -lm -L. -L../lib/ 
#LDFLAGS = -lgmp -lm -lmpfr -lgmpla -L.  -L../lib/
#LDFLAGS = -static -L./gmp-4.1.4/.libs/ -lgmp -lm  -L../lib/

SRCS = main.c convolution_multi_bignum.c convolution_multi_exact.c population.c
TARGET_DIR = ../target

UTILOBJS = ../target/util/factorial.o ../target/util/nck.o ../target/util/nckinit.o ../target/util/multichoose.o ../target/util/randi.o ../target/util/sum.o ../target/util/nnz.o ../target/util/nnzposcmp.o ../target/util/sortbynnzpos.o ../target/util/readmodel.o ../target/util/printmodel.o
GMPLAOBJS = ../target/gmpla/int_mat.o ../target/gmpla/int_matmatchrow.o ../target/gmpla/int_matprint.o ../target/gmpla/int_vec.o ../target/gmpla/int_vecprint.o ../target/gmpla/mpq_lubksb.o ../target/gmpla/mpq_ludcmp.o ../target/gmpla/mpq_matprint.o ../target/gmpla/mpq_matzeros.o ../target/gmpla/mpq_msp.o ../target/gmpla/mpq_mspget.o ../target/gmpla/mpq_mspprint.o ../target/gmpla/mpq_mspset_si.o ../target/gmpla/mpq_mspvecmul.o ../target/gmpla/mpq_vec.o ../target/gmpla/mpq_vecdup.o ../target/gmpla/mpq_vecprint.o

RM = /bin/rm -f

all: ca

ca: ${OBJS} util
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/ca ${LDFLAGS} 

util:  ${UTILOBJS}

# Include automatic dependency generation
include ../Makefile.deps

cleanbase: 
	${RM} ${OBJS}

clean: 
	-${RM} -rf $(TARGET_DIR)
	${RM} ${OBJS} ${DEPS} 
