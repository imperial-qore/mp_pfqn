# Multiprecision Product-Form Queueing Networks Solvers

The mp_pfqn library offers fast C solvers for product-form queueing networks using the GNU Multiple Precision Arithmetic Library (GMP). Exact arithmetic in GMP allows to compute normalizing constants for the equilibrium state probabilities on large closed multiclass models. 

- Closed network solvers:
  - Mean Value Analysis (MVA) [1]
  - Convolution Algorithm (CA) [2]
  - Recursion by Chain Algorithm (RECAL) [3]
  - Method of Moments (MoM) [4]
  - Class-Oriented Method of Moments (CoMoM) [5]
  - Generalized Method of Moments (gmom, divide-and-conquer, b=1) [9]
- Load-dependent solvers:
  - Generalized Load-Dependent (GLD) normalizing constant [6]
  - CoMoM-LD for load-dependent repairman models
- Mixed open/closed network solvers:
  - MVA-MX for mixed networks [7]
  - MVA-LD-MX for mixed networks with load-dependent stations [8]

## Quick Start
Run the following commands from the project root folder:
```bash
# Install all dependencies automatically (requires sudo)
sudo make install-deps

# Build the project
source ./setup-env.sh
make

# Run a solver
./bin/mva models/02_bottleneck_study.qn
./bin/ca models/02_bottleneck_study.qn
./bin/recal models/02_bottleneck_study.qn
./bin/mom models/02_bottleneck_study.qn
./bin/comom models/02_bottleneck_study.qn
```
If you wish to recompile after changes, run 
```
make clean
make
```

## Usage

### Model File Format

Models are specified in `.qn` text files with the following format:

```
R                   # Number of job classes
N1 N2 ... NR        # Population per class
Z1 Z2 ... ZR        # Think times per class
M                   # Number of queueing stations
m1 L11 L12 ... L1R  # Queue 1: multiplicity, demands at station 1 per class
m2 L21 L22 ... L2R  # Queue 2: multiplicity, demands at station 2 per class
...
mM LM1 LM2 ... LMR  # Queue M: multiplicity, demands at station M per class
```
The demand Lir=Sir*Vir is the product of the mean service time Sir and mean number of visits Vir, both for class-r jobs at station i. See the routing2visits utility to convert a routing probability matrix into the corresponding Vir values.

Optional keyword sections can follow the station lines:
```
LAMBDA               # Open-class arrival rates (use N=-1 for open classes)
l1 l2 ... lR         # Arrival rates per class (rational values allowed)
MU                   # Load-dependent service rates
mu_11 mu_12 ... mu_1Nt   # Station 1: rates for 1..Nt jobs
...
mu_M1 mu_M2 ... mu_MNt   # Station M: rates for 1..Nt jobs
```

*IMPORTANT*: All values must be passed as integers or rationals (e.g., `1/10`).

Examples are available under the models/ folder.

### Running Solvers

```bash
# Generate a random model
./bin/rndmodel M R N seed > model.qn
# Where: M = number of queues, R = number of classes, 
#        N = total population, seed = random seed

# Run different solvers (closed networks)
./bin/mva model.qn      # Mean Value Analysis
./bin/ca model.qn       # Convolution Algorithm
./bin/recal model.qn    # Recursion by Chain Algorithm
./bin/mom model.qn      # Method of Moments
./bin/comom model.qn    # Class-Oriented Method of Moments
./bin/gmom model.qn     # Generalized (divide-and-conquer) Method of Moments

# Load-dependent solvers (closed networks with MU section)
./bin/gld model.qn      # GLD normalizing constant
./bin/comomld model.qn  # CoMoM-LD for repairman models

# Mixed open/closed network solvers (models with LAMBDA section)
./bin/mvamx model.qn    # Mixed MVA
./bin/mvaldmx model.qn  # Mixed load-dependent MVA (requires MU section)
```

### Solver Options

All solvers support various command-line options to control their output format and behavior.

#### MVA Solver Options

```bash
./bin/mva [options] models/02_bottleneck_study.qn
```

| Option | Long Form | Description |
|--------|-----------|-------------|
| `-l` | `--log` | Print only log of normalizing constant as double |
| `-e` | `--ex` | Print exact normalizing constant numerator and denominator |
| `-g` | `--nc` | Print normalizing constant as double |
| `-t` | `--tput` | Print only throughputs, one per row |
| `-q` | `--qlen` | Print only queue lengths, one per row |
| `-d` | `--exact` | Print all performance metrics in full exact precision (integer or rational) |
| `-h` | `--help` | Print help message |


#### CA Solver Options

```bash
./bin/ca [options] models/02_bottleneck_study.qn
```

CA supports the same command-line options as MVA (see table above).


#### RECAL Solver Options

```bash
./bin/recal [options] models/02_bottleneck_study.qn
```

RECAL supports the same command-line options as MVA (see table above).


#### MoM Solver Options

```bash
./bin/mom [options] models/02_bottleneck_study.qn
```

| Option | Long Form | Description |
|--------|-----------|-------------|
| `-l` | `--log` | Print only log of normalizing constant as double |
| `-e` | `--ex` | Print exact normalizing constant numerator and denominator |
| `-g` | `--nc` | Print normalizing constant as double |
| `-t` | `--tput` | Print only throughputs, one per row |
| `-q` | `--qlen` | Print only queue lengths, one per row |
| `-d` | `--exact` | Print all performance metrics in full exact precision (integer or rational) |
| `-p digit` | | Apply perturbation at the specified digit (e.g., `-p 5`) |
| `-s seed` | `--seed` | Random seed for perturbation (default: 23000) |
| `-h` | `--help` | Print help message |

#### CoMoM Solver Options

```bash
./bin/comom [options] models/02_bottleneck_study.qn
```

CoMoM supports the same command-line options as MoM (see table above). The Class-Oriented Method of Moments provides an efficient algorithm specifically designed for multiclass models with large population sizes.

#### gmom Solver Options

```bash
./bin/gmom [ -e | -g | -l | -t | -q | --validate ] models/02_bottleneck.qn
```

`gmom` implements the generalized (divide-and-conquer, b=1) Method of
Moments: the model on the prefix `{1..m}` of queues is built from the model
on `{1..m-1}` plus one queue-removal (subtractive-convolution) branch, and
the overdetermined per-level system is solved by exact normal equations.
It supports `-e` (exact `G(N)` num/den), `-g` (`G(N)` double), `-l`
(`log G`), `-t` (throughputs), `-q` (queue lengths), and `--validate`
(check the whole moment basis against exact convolution). It requires at
least two classes and closed, distinct single-server stations; open/mixed,
load-dependent (`MU`), and multiserver (`mi>1`) models are rejected. `G`,
`X`, and `Q` match `ca` bit-for-bit on the supported models.

#### Examples

```bash
# Get only throughputs
./bin/mom -t models/02_bottleneck_study.qn

# Get only queue lengths  
./bin/mom -q models/02_bottleneck_study.qn

# Get exact normalizing constant (numerator and denominator)
./bin/mom -e models/02_bottleneck_study.qn

# Get all performance metrics in exact rational form
./bin/mom -d models/02_bottleneck_study.qn
./bin/mva -d models/02_bottleneck_study.qn

# Apply perturbation for approximate solution (MoM/CoMoM)
./bin/mom -p 5 models/02_bottleneck_study.qn
./bin/comom -p 5 models/02_bottleneck_study.qn

# Apply perturbation with specific random seed (MoM/CoMoM)
./bin/mom -p 5 -s 12345 models/02_bottleneck_study.qn
./bin/comom -p 5 -s 12345 models/02_bottleneck_study.qn
```
**Notes**: 
- The `-d`/`--exact` option is available in all solvers and prints performance metrics (throughputs, queue lengths, normalizing constants) in full exact precision as rational numbers (e.g., `1/3` instead of `0.333333`). This is useful for verification and when exact values are required.
- The perturbation option (`-p`) is available in both the MoM and CoMoM solvers and is useful when the exact solution fails due to singular systems. It introduces small randomized numerical perturbations to enable approximate solutions. When perturbation is applied, model parameters are displayed as original integer values with separate perturbation annotations (e.g., `50 eps=1.0e-05`). The `-s` option allows you to specify a random seed for reproducible perturbations.

## References

[1]: Reiser & Lavenberg (1980), *Mean-Value Analysis of Closed Multichain Queuing Networks,* Journal of the ACM 27(2).

[2]: Reiser, M. and Kobayashi, H. (1974), *Queuing Networks with Multiple Closed Chains: Theory and Computational Algorithms,* IBM Research Report RC-4919, July, 1974.

[3]: Conway, A. E. and Georganas, N. D. (1986), *RECAL—a new efficient algorithm for the exact analysis of multiple-chain closed queuing networks,* J. ACM 33, 4 (Oct. 1986), 768–791.

[4]: Casale (2006), *An efficient algorithm for the exact analysis of multiclass queueing networks with large population sizes,* Proc. of ACM SIGMETRICS 2006.

[5]: Casale, G. (2009), *CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models,* IEEE Trans. Softw. Eng. 35, 2 (March 2009), 162–177.

## License

This software is provided as-is for academic and research purposes and released as BSD-3.
