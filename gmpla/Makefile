CC = gcc
LD = gcc
CXX = g++
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../gmpla/ -I../util/ -Wall

# LinBox C++ flags (only used for linbox_solve.cpp)
LINBOX_CXXFLAGS = -std=c++11 -g -O2 -I../deps/local/include -I. -I../include/ -Wall -Wno-unused-parameter -Wno-deprecated-copy

SRCS = int_vecadd.c int_vecsubsum.c int_vecdiff.c int_mat.c int_matmatchrow.c int_matprint.c int_vec.c int_vecprint.c mpq_lubksb.c mpq_ludcmp.c mpq_matprint.c mpq_matzeros.c mpq_msp.c mpq_mspget.c mpq_mspprint.c mpq_mspset_si.c mpq_mspset_z.c mpq_mspvecmul.c mpq_vec.c mpq_vecdup.c mpq_vecprint.c mpq_matset_si.c mpq_matvecmul.c mpq_mattransmul.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

# LinBox C++ object (conditionally built)
LINBOX_OBJ = linbox_solve.o

RM = /bin/rm -f

gmpla: ${OBJS} ${LINBOX_OBJ}

# C++ compilation rule for LinBox wrapper
linbox_solve.o: linbox_solve.cpp linbox_solve.h
	$(CXX) $(LINBOX_CXXFLAGS) -c $< -o $@

# Include automatic dependency generation
include ../Makefile.deps

clean:
	${RM} ${OBJS} ${LINBOX_OBJ} ${DEPS}
