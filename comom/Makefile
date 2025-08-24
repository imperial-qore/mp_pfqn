CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 
#LDFLAGS = -lgmp -lm -L. -L../lib/ 
LDFLAGS = -L../deps/local/lib -lgmp -lmpfr -lm -L. -L../lib/
#LDFLAGS = -static -L./gmp-4.1.4/.libs/ -lgmp -lm  -L../lib/

SRCS = main.c setupls.c hash.c perfindices.c
TARGET_DIR = ../target

UTILOBJS = ../target/util/factorial.o ../target/util/nck.o ../target/util/nckinit.o ../target/util/multichoose.o ../target/util/randi.o ../target/util/sum.o ../target/util/nnz.o ../target/util/nnzposcmp.o ../target/util/sortbynnzpos.o ../target/util/readmodel.o ../target/util/printmodel.o
GMPLAOBJS = ../target/gmpla/int_mat.o ../target/gmpla/int_matmatchrow.o ../target/gmpla/int_matprint.o ../target/gmpla/int_vec.o ../target/gmpla/int_vecprint.o ../target/gmpla/mpq_lubksb.o ../target/gmpla/mpq_ludcmp.o ../target/gmpla/mpq_matprint.o ../target/gmpla/mpq_matzeros.o ../target/gmpla/mpq_msp.o ../target/gmpla/mpq_mspget.o ../target/gmpla/mpq_mspprint.o ../target/gmpla/mpq_mspset_si.o ../target/gmpla/mpq_mspvecmul.o ../target/gmpla/mpq_vec.o ../target/gmpla/mpq_vecdup.o ../target/gmpla/mpq_vecprint.o ../target/gmpla/int_vecadd.o ../target/gmpla/int_vecdiff.o ../target/gmpla/int_vecsubsum.o ../target/gmpla/mpq_matset_si.o 

RM = /bin/rm -f

all: comom

comom: ${OBJS} util
	${LD} ${OBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/comom ${LDFLAGS} 

util:  ${UTILOBJS}

# Include automatic dependency generation
include ../Makefile.deps

cleanbase: 
	${RM} ${OBJS}

clean: 
	-${RM} -rf $(TARGET_DIR)
	${RM} ${OBJS} ${DEPS} 
