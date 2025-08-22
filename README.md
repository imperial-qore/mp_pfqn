# Product-Form Queueing Networks (PFQN) Solvers

A C library for solving product-form queueing networks using exact multiprecision arithmetic. The latter allows to compute normalizing constants of the state probabilities in large models. 

## Features

- **Multiple Solving Algorithms**:
  - Mean Value Analysis (MVA)
  - Method of Moments (MoM)

## Prerequisites

- C/C++ compiler (gcc/g++ recommended)
- GNU Make

## Quick Start

```bash
# Install all dependencies automatically (requires sudo)
sudo make install-deps

# Build the project
make

# Run a solver
./bin/mva models/02_bottleneck_study.qn
./bin/mom models/02_bottleneck_study.qn
```

## Usage

### Model File Format

Models are specified in `.qn` text files with the following format:

```
R                    # Number of classes
N1 N2 ... NR        # Population per class
Z1 Z2 ... ZR        # Think times per class
M                    # Number of queues
m1 L11 L12 ... L1R  # Queue 1: multiplicity, demands per class
m2 L21 L22 ... L2R  # Queue 2: multiplicity, demands per class
...
mM LM1 LM2 ... LMR  # Queue M: multiplicity, demands per class
```
*IMPORTANT*: All values must be passed as integers. 

### Running Solvers

```bash
# Generate a random model
./bin/rndmodel M R N seed > model.qn
# Where: M = number of queues, R = number of classes, 
#        N = total population, seed = random seed

# Run different solvers
./bin/mva model.qn      # Mean Value Analysis
./bin/mom model.qn      # Method of Moments
```

## References

[1]: Reiser & Lavenberg (1980), *Mean-Value Analysis of Closed Multichain Queuing Networks,* Journal of the ACM 27(2).

[2]: Casale (2006), *An efficient algorithm for the exact analysis of multiclass queueing networks with large population sizes,* Proc. of ACM SIGMETRICS 2006.

## License

This software is provided as-is for academic and research purposes and released as BSD-3.
