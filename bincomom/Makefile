CC = gcc
LD = gcc
CFLAGS = -g -O4 -fopenmp -I../deps/local/include -I. -I../comom/ -I../include/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -fopenmp -L../deps/local/lib -lgmp -lmpfr -lm -L. -L../lib/

SRCS = main.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

# bin/comom itself: its main.c compiled here with main renamed (bincomom's
# main sets the COMOM_ADVANCE hook and calls it), its other objects and the
# comommod engine (tier 4 of the singular fallback) taken from ../comom
COMOMOBJS = ../comom/setupls.o ../comom/hash.o ../comom/mpq_ludcmp_progress.o ../comom/btf_decompose.o ../comom/btf_solve.o ../comom/ccengine.o

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_ludcmp.o ../gmpla/mpq_matprint.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_msp.o ../gmpla/mpq_mspget.o ../gmpla/mpq_mspprint.o ../gmpla/mpq_mspset_si.o ../gmpla/mpq_mspset_z.o ../gmpla/mpq_mspvecmul.o ../gmpla/mpq_vec.o ../gmpla/mpq_vecdup.o ../gmpla/mpq_vecprint.o ../gmpla/int_vecadd.o ../gmpla/int_vecdiff.o ../gmpla/int_vecsubsum.o ../gmpla/mpq_matset_si.o

RM = /bin/rm -f

all: bincomom

bincomom: ${OBJS} comom_main.o comomobjs util
# -fopenmp on the link line itself: the top-level Makefile passes LDFLAGS on
# the command line, which replaces the assignment above and with it -fopenmp
	${LD} ${OBJS} comom_main.o ${COMOMOBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/bincomom ${LDFLAGS} -fopenmp

comom_main.o: ../comom/main.c ../comom/comom.h
	${CC} ${CFLAGS} -Dmain=comom_cli_main -c $< -o $@

.PHONY: comomobjs
comomobjs:
	$(MAKE) -C ../comom $(notdir $(COMOMOBJS))

util:  ${UTILOBJS}

include ../Makefile.deps

cleanbase:
	${RM} ${OBJS} comom_main.o

clean:
	${RM} ${OBJS} comom_main.o ${DEPS}
