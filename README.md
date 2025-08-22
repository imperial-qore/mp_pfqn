# Product-Form Queueing Networks (PFQN) Solvers

A comprehensive C library for solving product-form queueing networks using exact multiprecision arithmetic for computing normalizing constants. 

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
```

## Build Commands

```bash
# Build all programs
make clean; make
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

### Example

Create a simple model file `example.qn`:
```
2          # 2 job classes
10 5       # 10 jobs of class 1, 5 jobs of class 2
0 0        # No think times
3          # 3 queues
1 10 5     # Queue 1 (multiplicity 1, class-1 demand = 10, class-2 demand = 5)
1 5 9      # Queue 2
1 11 7     # Queue 3
```

Run the MVA solver:
```bash
./bin/mom example.qn
```

## Solvers & References

### MVA (Mean Value Analysis)
Classic iterative algorithm for product-form networks. Efficient for small to medium-sized models.  
*Reference: Reiser & Lavenberg (1980), "Mean-Value Analysis of Closed Multichain Queuing Networks," Journal of the ACM 27(2).*
```bash
./bin/mva models/02_bottleneck_study.qn
```

### MoM (Method of Moments)
Normalizing constant-based approach using moment generating functions. Suitable for models with moderate populations.  
*Reference: Casale (2006), "An efficient algorithm for the exact analysis of multiclass queueing networks with large population sizes," Proc. of ACM SIGMETRICS 2006.*
```bash
./bin/mom models/02_bottleneck_study.qn
```

## Performance Considerations

- For small models (< 50 jobs), MVA is typically fastest
- For larger populations, MoM may be more efficient
- The MoM algorithm may fail on particular combinations of demands due to the underpinning systems of linear equations becoming singular. Use -p to introduce perturbations on the demand and solve an approximate model. 

## License

This software is provided as-is for academic and research purposes and released as BSD-3.
