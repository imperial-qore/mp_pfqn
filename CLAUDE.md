# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**mp_pfqn** is a C library for exact analysis of Product-Form Queueing Networks using GMP (GNU Multiple Precision) arithmetic. It computes normalizing constants, throughputs, and queue lengths for closed multiclass models. All arithmetic is done in exact rational precision via `mpq_t`/`mpz_t` types.

## Build & Run

```bash
# First time: install dependencies (GMP, MPFR)
sudo make install-deps

# Build all solvers
source ./setup-env.sh
make

# Rebuild after changes
make clean && make
```

Binaries are placed in `bin/`: `mva`, `ca`, `recal`, `mom`, `comom`, `procomom`, `rndmodel`, `routing2visits`.

## Testing

```bash
# Run full regression suite (compares output parity across all solvers)
./run_unit_tests.sh

# Run one model with all options
./run_unit_tests.sh models/03_think.qn

# Run one model with a specific option
./run_unit_tests.sh models/03_think.qn -- -t

# Run multiple specific options
./run_unit_tests.sh models/02_bottleneck.qn -- -e -l
```

The test script runs all 5 solvers (mva, ca, recal, mom, comom) with options `-e`, `-l`, `-t`, `-q` on each `.qn` model file, then diffs outputs for parity. Timeout is 120s per test. Errors are logged to `test_errors_*.log`.

## Architecture

### Solver Modules (each has its own `main.c` and `Makefile`)

| Directory | Algorithm | Notes |
|-----------|-----------|-------|
| `mva/` | Mean Value Analysis | Classical iterative approach |
| `ca/` | Convolution Algorithm | With multi-server station expansion |
| `recal/` | Recursion by Chain (RECAL) | Efficient chain-based recursion |
| `mom/` | Method of Moments (MoM) | Exact linear system solve; optional LinBox C++ backend |
| `comom/` | Class-Oriented MoM (CoMoM) | BTF decomposition for sparse systems |
| `procomom/` | ProCoMoM | Marginal queue-length probabilities |
| `clw/` | Choudhury-Leung-Whitt (CLW) | Generating-function inversion; single-server + IS. Floating-point (G/logG as doubles, not exact) |
| `clwld/` | CLW with load-dependent stations | Bertozzi-McKenna per-center transforms; multiserver and load-dependent. Floating-point output |

### Shared Libraries

- **`gmpla/`** — GMP linear algebra: dense matrices (`mpq_mat_t`), sparse matrices (`mpq_msp_t`), LU decomposition (`mpq_ludcmp`/`mpq_lubksb`), LinBox interface
- **`util/`** — Model I/O (`readmodel`), population enumeration (`nextpop`, `popindex`), combinatorics (`nck`)

### Core Data Structure (`util/util.h`)

```c
typedef struct {
    int M;           // number of queues
    int R;           // number of classes
    int *N;          // job populations per class
    mpz_t *Z;        // think times (GMP integers)
    mpz_t **L;       // service demands (GMP integers)
    int *mi;         // queue multiplicities (server count)
} qnmodel;
```

### Model File Format (`.qn`)

All values must be integers. Rationals are encoded via GMP's `num/den` syntax.

```
R                    # number of job classes
N1 N2 ... NR         # population per class
Z1 Z2 ... ZR         # think times per class
M                    # number of stations
m1 L11 L12 ... L1R   # station 1: multiplicity, demands per class
...
```

12 test models in `models/` (01_single through 12_expanded).

### Solver CLI Options

All solvers share: `-l` (log NC), `-e` (exact NC num/den), `-g` (NC as double), `-t` (throughputs), `-q` (queue lengths), `-d` (all metrics exact rational).

MoM/CoMoM additionally support: `-p digit` (perturbation), `-s seed` (perturbation seed).

## Build System Details

- Top-level `Makefile` orchestrates builds of all subdirectories sequentially
- Each subdirectory has its own `Makefile` that compiles `.c` files and links against shared `util/` and `gmpla/` objects
- `mom/` links with `g++` (for LinBox C++ integration); all others use `gcc`
- Compiler flags: `-g -Wall -O4` plus GMP/MPFR include/link paths
- Dependencies can live in `deps/local/` (manual build) or system paths

## Commit Convention

Use conventional commits: `feat:`, `fix:`, `chore:`, etc.
