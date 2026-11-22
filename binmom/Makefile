CC = gcc
LD = gcc
CFLAGS = -g -O4 -fopenmp -I../deps/local/include -I. -I../mom/ -I../include/ -I../gmpla/ -I../util/ -Wall
LDFLAGS = -fopenmp -L../deps/local/lib -lgmp -lm -L. -L../lib/

SRCS = main.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

# mom's exact setup, block solve, descent and reports, compiled here so that
# bin/mom's build is untouched (no LinBox: binmom never sets USE_LINBOX)
MOMSRCS = setupls.c blocksolve.c mdecrease.c perfindices.c mpq_ludcmp_progress.c
MOMOBJS = $(addprefix mom_,$(MOMSRCS:.c=.o))

UTILOBJS = ../util/factorial.o ../util/nck.o ../util/nckinit.o ../util/multichoose.o ../util/randi.o ../util/sum.o ../util/nnz.o ../util/nnzposcmp.o ../util/sortbynnzpos.o ../util/readmodel.o ../util/printmodel.o
GMPLAOBJS = ../gmpla/int_mat.o ../gmpla/int_matmatchrow.o ../gmpla/int_matprint.o ../gmpla/int_vec.o ../gmpla/int_vecprint.o ../gmpla/mpq_lubksb.o ../gmpla/mpq_ludcmp.o ../gmpla/mpq_matprint.o ../gmpla/mpq_matzeros.o ../gmpla/mpq_msp.o ../gmpla/mpq_mspget.o ../gmpla/mpq_mspprint.o ../gmpla/mpq_mspset_si.o ../gmpla/mpq_mspset_z.o ../gmpla/mpq_mspvecmul.o ../gmpla/mpq_vec.o ../gmpla/mpq_vecdup.o ../gmpla/mpq_vecprint.o

RM = /bin/rm -f

all: binmom

binmom: ${OBJS} ${MOMOBJS} util
# -fopenmp on the link line itself: the top-level Makefile passes LDFLAGS on
# the command line, which replaces the assignment above and with it -fopenmp
	${LD} ${OBJS} ${MOMOBJS} ${GMPLAOBJS} ${UTILOBJS} -o ../bin/binmom ${LDFLAGS} -fopenmp

mom_%.o: ../mom/%.c
	${CC} ${CFLAGS} -c $< -o $@

util:  ${UTILOBJS}

include ../Makefile.deps

cleanbase:
	${RM} ${OBJS} ${MOMOBJS}

clean:
	${RM} ${OBJS} ${MOMOBJS} ${DEPS}
