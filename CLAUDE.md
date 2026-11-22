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

Binaries are placed in `bin/`: `mva`, `ca`, `recal`, `mom`, `mommod`, `comom`, `comommod`, `procomom`, `rndmodel`, `routing2visits`.

**`bool` in headers**: `util/util.h` includes `<stdbool.h>` unconditionally. It must
not go back to defining `bool` as `unsigned int`: sources differ in whether they
include `<stdbool.h>`, so that made the same global a 1-byte `_Bool` in one object
and a 4-byte `unsigned int` in another, and `mom` silently returned wrong results in
a few percent of runs. See `mommod/README.md`.

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
| `mom/` | Method of Moments (MoM) | Exact linear system solve; optional LinBox C++ backend. Singular models are solved exactly by the tiered fallback of `util/singular.h` instead of the old automatic perturbation; tier 4 runs on mommod's modular engine (`mommod/engine.c`, linked as `mom/mmengine.o` with every symbol but the `mmengine.h` API made local), `-j` threads |
| `comom/` | Class-Oriented MoM (CoMoM) | BTF decomposition for sparse systems. Singular models are solved exactly by the tiered fallback of `util/singular.h` (derived models re-solved through `comom_main` in capture mode); tier 4 runs on comommod's modular engine (`comommod/engine.c`, linked as `comom/ccengine.o` with every symbol but the `ccengine.h` API made local), `-j` threads; `-X` fails instead, `-p` keeps the approximate perturbation |
| `gmom/` | Generalized MoM (divide-and-conquer, b=1) | Prefix-chain queue-removal branch (one queue per step); overdetermined per-level system solved by exact normal equations; lower replica levels carried through the recursion give G(N)/X/Q; distinct single-server closed queues only. `-p`/`-s` perturbation with analytic descaling, `-b` two-sided (+eps/-eps) bracket, `--force-perturb`, `-X`. See `gmom/README.md` |
| `safe_comom/` | CoMoM exact on singular models | Hybrid MVA/CoMoM (TSE'09 Appendix A): tiered orchestrator - plain CoMoM, then class-permuted CoMoM, then exact convolution (`ca`) for genuine loading degeneracies. Every tier exact. See `safe_comom/README.md` |
| `safe_mom/` | MoM exact on singular models | Same tiered construction on the MoM basis: plain `mom -X`, then class-permuted `mom -X`, then exact convolution (`ca`). Adds `-X`/`--no-perturb` to `mom` (fail instead of auto-perturbing). Tier driver shared with `safe_comom` (`util/safe_tiers.c`). See `safe_mom/README.md` |
| `safe_gmom/` | Generalized MoM exact on singular models | Same tiers on gmom's b=1 basis: `gmom -X`, class-permuted `gmom -X`, then `ca`. Adds `-X`/`--no-perturb` to `gmom`. Tier 2 never fires (gmom picks the class order itself via its exact b=1 oracle); Tier 3 also absorbs gmom's domain limits (R=1, multiserver, load-dependent). Open/mixed reported as out of scope. See `safe_gmom/README.md` |
| `procomom/` | ProCoMoM | Marginal queue-length probabilities on the CoMoM basis. Singular systems solved exactly by `util/singular_marg.h` (class orders, merged stations/classes, interpolation); `-p` keeps the old approximate perturbation, `-X` fails instead. Three defects fixed after parity testing against `promom`/`ca`: missing `mi` weights in the PC rows (was wrong at every station of any model with a replicated station); no basis column for the reference station's own replica (stride now `M+1` only when that station is replicated -- widening it unconditionally makes `A^T A` singular); and duplicated output rows on the perturbation retry. `-P` now matches `promom` bit-for-bit. See `procomom/README.md` `--progress` reports per-step timing on stderr; two cost regimes, large basis (AtA products, ~75% of a step, recomputed each step though almost nothing changes) and post-perturbation coefficient growth (LU 0.01s -> 143s on `05_sparse`, impractical -- use `camarg`). |
| `promom/` | ProMoM | Marginal queue-length probabilities on the MoM basis (replica combinations x class decrements) instead of the CoMoM one. Same `A q(n) = B q_prev(n) + n[DC q(n-1) + DD q_prev(n-1)]` shape as procomom; the population constraint holds the reference copy out of the station sum. Singular systems solved exactly by `util/singular_marg.h` (`-X` fails instead, used by `safe_promom`). `-q` matches exact `ca` to 5e-16; `-P` is bit-for-bit `procomom` on single-server models and CORRECT where procomom is wrong (any `mi>1`). See `promom/README.md` |
| `camarg/` | Convolution marginals | Exact marginal queue-length distribution per station by Buzen convolution over the population lattice, `q_k(N,n) = sum_{|j|=n} (n!/prod j_r!) prod L_kr^{j_r} G_-k(N-j)`. No coefficient matrix, so never singular: the oracle for `promom`/`procomom` on degenerate models and `safe_promom`'s Tier 3. Cost is exponential in R. See `camarg/README.md` |
| `safe_promom/` | ProMoM exact on class-order singularities | Tier wrapper over `util/safe_tiers`: `promom -X`, then class-permuted `promom -X`. Tier 3 is `camarg` (exact convolution marginals), not `ca`, which returns means rather than distributions. Tier 2 recovers `08_multiclass` and `11_swapped`. See `safe_promom/README.md` |
| `binmom/` | MoM by binary splitting | Same recursion and output as `mom` (byte-identical), but within each class the steps `n g(n) = (nP+Q) g(n-1)` (blocks fixed per class) are multiplied together by binary splitting over the integers: about `d^3 N polylog N` instead of `d^2 N^2 log N`. Splits a class only when `N_r >= 2 d^2` (`--cross`), else takes `mom`'s steps; `-j` threads the products. Up to 500x faster than every other solver at small basis and large population; matches `mom` elsewhere; `mommod` stays faster at moderate basis. Prototype: refuses mixed zero/non-zero think times and singular systems (use `mom`). See `binmom/README.md` |
| `bincomom/` | CoMoM by binary splitting | `bin/comom` itself (its `main.c` compiled with `main` renamed, comom's objects and `ccengine.o` linked) plus one hook, `COMOM_ADVANCE`: once a class system is set up it moves the basis from `n_r = 0` to `N_r - (M+1)` by binary splitting (`util/binsplit.h`); the last `M+1` steps stay with comom. The step matrices are read off `comom_step()` (comom's step, moved unchanged into a function) at `n = 1, 2`, so every special path of comom is inside them; a guard checks them against a real step. Byte-identical to `comom`, every flag and the singular fallback included, no refusals of its own. Splits when `N_r - M - 1 >= 2 d^2` (`--cross`). 330x faster than comom at small basis and large population, the same elsewhere. See `bincomom/README.md` |
| `mommod/` | Multi-modular MoM | Whole recursion replayed in `Z/p` for many word-size primes, CRT plus rational reconstruction at the end, OpenMP over primes (`-j`). `-B` halves the prime count via denominator-bounded reconstruction, verified against witness primes held out of the CRT (`-W`); off by default since it is probabilistic. Singular models are solved exactly by a tiered fallback: other class orders, merging identical stations or proportional classes, then interpolation of the exact limit of a perturbation (tiers 1-3 shared in `util/singular.h`, tier 4 per prime). See `mommod/README.md` |
| `comommod/` | Multi-modular CoMoM | The `comom` class-oriented recursion replayed in `Z/p` for many word-size primes, CRT plus rational reconstruction at the end, OpenMP over primes (`-j`). Same flags as `mommod`, except that `-B`'s denominator bound defaults to a third of the modulus: CoMoM denominators grow with the population instead of staying a few dozen bits wide. Singular models are solved exactly by the same tiered fallback as `mommod`. See `comommod/README.md` |

### Shared Libraries

- **`gmpla/`** — GMP linear algebra: dense matrices (`mpq_mat_t`), sparse matrices (`mpq_msp_t`), LU decomposition (`mpq_ludcmp`/`mpq_lubksb`), LinBox interface
- **`zpla/`** — the same API over `Z/p` in machine words (`zp_*`), Montgomery arithmetic, CRT and rational reconstruction
- **`util/`** — Model I/O (`readmodel`), population enumeration (`nextpop`, `popindex`), combinatorics (`nck`),
  the shared b=1 level builder and singularity oracle (`setup1`, `pfqn_recursion_singular`,
  `pfqn_nonsingular_recclass` in `pfqn_sing.{c,h}`), the tiered exact-wrapper driver
  (`safe_tiers.{c,h}`, used by `safe_comom`, `safe_mom`, `safe_gmom` and `safe_promom`); header-only
  `lvtab.h` (replica-level tables of the full MoM basis, `mom_mixed_think`), `singular.h` (tiered exact
  fallback for singular MoM/CoMoM systems, shared by `mom`, `mommod`, `comom`, `comommod`) and
  `singular_marg.h` (the same tiers for per-station marginal distributions, interpolating the unnormalized
  marginals; used by `promom` and `procomom`), `binsplit.h` (exact products of a class's population steps by
  binary splitting, shared by `binmom` and `bincomom`) and `replicas.h` (the copies meaning of `mi` for the
  load-dependent solvers: `qn_replica_mu` folds `mi` copies into one load-dependent station, exactly;
  `qn_expand_replicas` writes one station per copy)

### Core Data Structure (`util/util.h`)

```c
typedef struct {
    int M;           // number of queues
    int R;           // number of classes
    int *N;          // job populations per class
    mpz_t *Z;        // think times (GMP integers)
    mpz_t **L;       // service demands (GMP integers)
    int *mi;         // multiplicity: number of identical copies of the station (see below)
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

**Multiplicity `mi`** means `mi` identical copies of the station, each a single-server queue with the
station's demands (and rates, with a MU section), every job visiting each copy; queue lengths are reported as
the total over the copies. Every solver reads it this way. A queue with `c` servers sharing one line is
written with a MU section, rates `1 2 ... c c c`. `gld`, `comomld` and `mvaldmx` used to read `mi` as `c`
servers (rates `min(n, mi)`), a different model, and disagreed with `ca` by 1% to 80% on models with `mi > 1`;
they now fold the copies into one load-dependent station of rates `c(t-1)/c(t)`, `c` the `mi`-fold convolution
of `1/(mu(1)...mu(t))` (`t/(t+mi-1)` for plain queues), or, in `mvaldmx`, expand them (`util/replicas.h`), and
match `ca` exactly.

12 test models in `models/` (01_single through 12_expanded).

### Solver CLI Options

All solvers share: `-l` (log NC), `-e` (exact NC num/den), `-g` (NC as double), `-t` (throughputs), `-q` (queue lengths), `-d` (all metrics exact rational).

MoM/CoMoM additionally support: `-p digit` (perturbation), `-s seed` (perturbation seed).

Letters keep one meaning across solvers: `-b` bounds (`mom`, `gmom`), `-P` marginal/state probabilities
(`promom`, `procomom`, `camarg`, `comomld`), `-q` mean queue lengths, `-W` witness primes (`mommod`,
`comommod`). Renamed to get there, the old letter now failing with a pointer to the new name: prime bit size
`-b` -> `--primebits` (`mommod`, `comommod`), back probabilities `-b` -> `--back` (`lcfsmva`), state probabilities
`-q` -> `-P` (`comomld`, which has no mean queue lengths).
`mom`, `comom`, `gmom`, `promom` and `procomom` also support `-X`/`--no-perturb`: fail on a singular system instead of running the exact fallback (or, for `gmom`, auto-perturbing).

## Build System Details

- Top-level `Makefile` orchestrates builds of all subdirectories sequentially
- Each subdirectory has its own `Makefile` that compiles `.c` files and links against shared `util/` and `gmpla/` objects
- `mom/` links with `g++` (for LinBox C++ integration); all others use `gcc`
- `mom/` also links `mmengine.o`, built by its Makefile from `mommod/` objects (`make -C ../mommod <objs>`) plus `zpla.o` via `ld -r` and `objcopy --keep-global-symbol` (only the `mommod/mmengine.h` API stays global, since both solvers define `setupls`, `blocksolve`, `mdecrease`, ...); plain C, no extra C++ compile
- `comom/` likewise links `ccengine.o` from `comommod/` objects (`engine.o globals.o image.o setupls.o combs.o btf_*.o`) plus `zpla.o`, keeping only the `comommod/ccengine.h` API global (both solvers define `setupls`, `hash`, `btf_*`, `qnm`, ...), and links with `-fopenmp` for the engine's threads
- Compiler flags: `-g -Wall -O4` plus GMP/MPFR include/link paths
- Dependencies can live in `deps/local/` (manual build) or system paths

## Commit Convention

Use conventional commits: `feat:`, `fix:`, `chore:`, etc.

## Development-only solvers

`momf` (MoM in MPFR fixed precision, with its `fpla/` library) lives on the `momf-dev` branch only, outside the production releases cut from `main` and `dev`. `momf-dev` follows `dev` by merging `dev` into it; the commit that removed `momf` from `dev` is recorded there as merged (`git merge -s ours`), so those merges keep `momf`.
