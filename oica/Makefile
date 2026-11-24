CC = gcc
CFLAGS = -g -O2 -I../deps/local/include -Wall
LDFLAGS = -L../deps/local/lib -lgmp -lm

all: oica

oica: main.c
	${CC} ${CFLAGS} main.c -o ../bin/oica ${LDFLAGS}

clean:
	/bin/rm -f ../bin/oica
