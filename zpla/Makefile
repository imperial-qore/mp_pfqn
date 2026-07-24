CC = gcc
LD = gcc
CFLAGS = -g -O4 -I../deps/local/include -I. -I../include/ -I../zpla/ -I../util/ -Wall

SRCS = zpla.c
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)

RM = /bin/rm -f

.PHONY: all
all: ${OBJS}

include ../Makefile.deps

clean:
	${RM} ${OBJS} ${DEPS}
