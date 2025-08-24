CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall 

SRCS = int_vecadd.c int_vecsubsum.c int_vecdiff.c int_mat.c int_matmatchrow.c int_matprint.c int_vec.c int_vecprint.c mpq_lubksb.c mpq_ludcmp.c mpq_matprint.c mpq_matzeros.c mpq_msp.c mpq_mspget.c mpq_mspprint.c mpq_mspset_si.c mpq_mspset_z.c mpq_mspvecmul.c mpq_vec.c mpq_vecdup.c mpq_vecprint.c mpq_matset_si.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

RM = /bin/rm -f

gmpla: ${OBJS} 

# Include automatic dependency generation
include ../Makefile.deps

clean: 
	${RM} ${OBJS} ${DEPS} 
