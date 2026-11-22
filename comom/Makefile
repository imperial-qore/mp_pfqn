CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 
#LDFLAGS = -lgmp -lm -L. -L../lib/ 
LDFLAGS = -fopenmp -L../deps/local/lib -lgmp -lmpfr -lm -L. -L../lib/
#LDFLAGS = -static -L./gmp-4.1.4/.libs/ -lgmp -lm  -L../lib/

SRCS = main.c setupls.c hash.c mpq_ludcmp_progress.c btf_decompose.c btf_solve.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_ludcmp.o ../gmpla/mpq_matprint.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_msp.o ../gmpla/mpq_mspget.o ../gmpla/mpq_mspprint.o ../gmpla/mpq_mspset_si.o ../gmpla/mpq_mspset_z.o ../gmpla/mpq_mspvecmul.o ../gmpla/mpq_vec.o ../gmpla/mpq_vecdup.o ../gmpla/mpq_vecprint.o ../gmpla/int_vecadd.o ../gmpla/int_vecdiff.o ../gmpla/int_vecsubsum.o ../gmpla/mpq_matset_si.o 

RM = /bin/rm -f

all: comom

comom: ${OBJS} util ccengine.o
# -fopenmp on the link line itself: the top-level Makefile passes LDFLAGS on
# the command line, which replaces the assignment above and with it -fopenmp
	${LD} ${OBJS} ccengine.o ${GMPLAOBJS} ${UTILOBJS} -o ../bin/comom ${LDFLAGS} -fopenmp

# The multi-modular engine of comommod, used for tier 4 of the singular-model
# fallback (as bin/mom uses mommod's).  comommod's sources reuse names that
# comom defines too (setupls, hash, btf_*, qnm, ...), so they are linked into
# one relocatable object whose symbols are all made local except the
# ccengine.h API.
CCOBJS = ../comommod/engine.o ../comommod/globals.o ../comommod/image.o ../comommod/setupls.o ../comommod/combs.o ../comommod/btf_decompose.o ../comommod/btf_solve.o
CCAPI = comommod_config comommod_stats comommod_solve_tiered comommod_interp

.PHONY: ccobjs
ccobjs:
	$(MAKE) -C ../comommod $(notdir $(CCOBJS))
	$(MAKE) -C ../zpla

ccengine.o: ccobjs
	ld -r -o ccengine_raw.o $(CCOBJS) ../zpla/zpla.o
	objcopy $(addprefix --keep-global-symbol=,$(CCAPI)) ccengine_raw.o ccengine.o
	${RM} ccengine_raw.o

util:  ${UTILOBJS}

# Include automatic dependency generation
include ../Makefile.deps

cleanbase: 
	${RM} ${OBJS}

clean: 
	${RM} ${OBJS} ${DEPS} ccengine.o 
