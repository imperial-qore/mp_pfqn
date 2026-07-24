CC	= gcc
PRJCFLAGS= -g -Wall -O4
LD	= gcc
LDFLAGS	= 
AR	= ar
ARFLAGS	=
RANLIB	= ranlib
RM	= rm
ECHO	= echo
SHELL	= /bin/sh

# Dependency configuration
DEPS_DIR = $(shell pwd)/deps/local

# Check for dependencies: manual install first, then system
ifneq ($(wildcard $(DEPS_DIR)),)
    DEP_CFLAGS = -I$(DEPS_DIR)/include
    DEP_LDFLAGS = -L$(DEPS_DIR)/lib -lgmp -lmpfr -lm
    DEP_SOURCE = manual
else
    # Fall back to system libraries
    DEP_CFLAGS = 
    DEP_LDFLAGS = -lgmp -lmpfr -lm
    DEP_SOURCE = system
endif

PRJCFLAGS += $(DEP_CFLAGS)
LDFLAGS += $(DEP_LDFLAGS)

#.SILENT:

DIRS	= util gmpla fpla zpla mva mvamx mom recal rndmodel ca comom routing2visits procomom gld comomld mvaldmx lcfsmva clw clwld momf mommod gmom safe_comom
EXE	= ./bin/mom ./bin/momf ./bin/mommod ./bin/gmom ./bin/comom ./bin/safe_comom ./bin/procomom
LIBS	= -L. -lsub -lsuba -lsubsub

all : check-deps $(EXE) 

check-deps:
	@if [ "$(DEP_SOURCE)" = "manual" ]; then \
		if [ ! -f "$(DEPS_DIR)/include/gmp.h" ]; then \
			echo "WARNING: GMP header not found in $(DEPS_DIR)/include/"; \
		fi; \
		if [ ! -f "$(DEPS_DIR)/include/mpfr.h" ]; then \
			echo "WARNING: MPFR header not found in $(DEPS_DIR)/include/"; \
		fi; \
	elif [ "$(DEP_SOURCE)" = "system" ]; then \
		if [ ! -f "/usr/include/gmp.h" ] && [ ! -f "/usr/local/include/gmp.h" ] && ! find /usr/include -name "gmp.h" -print -quit 2>/dev/null | grep -q gmp.h; then \
			echo ""; \
			echo "ERROR: GMP development headers not found!"; \
			echo "Please install libgmp-dev:"; \
			echo "  Ubuntu/Debian: sudo apt-get install libgmp-dev"; \
			echo "  RHEL/CentOS: sudo yum install gmp-devel"; \
			echo "  macOS: brew install gmp"; \
			echo ""; \
			echo "Or run 'make install-deps' to build from source"; \
			exit 1; \
		fi; \
		if [ ! -f "/usr/include/mpfr.h" ] && [ ! -f "/usr/local/include/mpfr.h" ] && ! find /usr/include -name "mpfr.h" -print -quit 2>/dev/null | grep -q mpfr.h; then \
			echo ""; \
			echo "ERROR: MPFR development headers not found!"; \
			echo "Please install libmpfr-dev:"; \
			echo "  Ubuntu/Debian: sudo apt-get install libmpfr-dev"; \
			echo "  RHEL/CentOS: sudo yum install mpfr-devel"; \
			echo "  macOS: brew install mpfr"; \
			echo ""; \
			echo "Or run 'make install-deps' to build from source"; \
			exit 1; \
		fi; \
	fi

$(EXE) : 
	@if [ ! -f $(EXE) ]; then make -s check-deps; fi
	@mkdir -p bin
	cd gmpla; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd mva; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd mvamx; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd mom; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd recal; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd util; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd rndmodel; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd ca; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd comom; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd routing2visits; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd procomom; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd gld; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd comomld; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd mvaldmx; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd lcfsmva; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd clw; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	cd clwld; $(MAKE) $(MFLAGS) PRJCFLAGS="$(PRJCFLAGS)" LDFLAGS="$(LDFLAGS)"; cd ..
	@echo "Cleaning up .o files..."
	@find . -name "*.o" -type f -not -path "./deps/*" -delete

install-deps:
	@echo "Installing dependencies from source..."
	bash install-deps.sh

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force :
	true

.PHONY: all check-deps install-deps clean force

