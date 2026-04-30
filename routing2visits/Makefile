CC = gcc
CFLAGS = $(PRJCFLAGS) -std=c99
LIBS = -lgmp

SOURCES = routing2visits.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = routing2visits
LIBRARY = librouting2visits.a

.PHONY: all clean library

# routing2visits.c currently contains only library functions; the CLI
# driver (main()) was never committed, so the executable target is
# disabled and only the static library is built.
all: library

library: $(LIBRARY)

$(LIBRARY): routing2visits.o
	ar rcs $@ $^

$(EXECUTABLE): routing2visits.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)
	cp $(EXECUTABLE) ../bin/

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

routing2visits.o: routing2visits.c routing2visits.h

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(LIBRARY)

install: library
	@echo "To install, copy librouting2visits.a and routing2visits.h to appropriate directories"
	@echo "Example: cp $(LIBRARY) /usr/local/lib/"
	@echo "         cp routing2visits.h /usr/local/include/"

test: $(EXECUTABLE)
	@echo "To test the utility, create input files and run:"
	@echo "./$(EXECUTABLE) <num_stations> <num_classes> <routing_prefix> <service_times_file>"
	@echo "See example files in examples/ directory"

help:
	@echo "Available targets:"
	@echo "  all          - Build library and routing2visits utility (default)"
	@echo "  library      - Build static library librouting2visits.a"
	@echo "  routing2visits - Build routing2visits utility"
	@echo "  test         - Show test instructions"
	@echo "  clean        - Remove all built files"
	@echo "  install      - Show installation instructions"
	@echo "  help         - Show this help message"