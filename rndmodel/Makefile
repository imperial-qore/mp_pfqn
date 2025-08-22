CFLAGS = -g -Wall
CC = gcc
LD = gcc 
LDFLAGS = -lm
OBJS = rndmodel.o 
PROG = rndmodel 
RM = /bin/rm -f

all: $(PROG) 

$(PROG): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o ../bin/$(PROG)

clean: 
	$(RM) $(OBJS)
