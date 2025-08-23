# Product-Form Queueing Networks (PFQN) Solvers

The mp_pfqn library offers fast C solvers for product-form queueing networks using the GNU Multiple Precision Arithmetic Library (GMP). Exact arithmetic in GMP allows to compute normalizing constants for the equilibrium state probabilities on large closed multiclass models. 

- Current solvers:
  - Mean Value Analysis (MVA) [1]
  - Method of Moments (MoM) [2]

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
./bin/mom models/02_bottleneck_study.qn
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
*IMPORTANT*: All values must be passed as integers. 

Examples are available under the models/ folders.

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

### Solver Options

Both solvers support various command-line options to control their output format and behavior.

#### MVA Solver Options

```bash
./bin/mva [options] model.qn
```

| Option | Long Form | Description |
|--------|-----------|-------------|
| `-v` | `--verbose` | Print exact ratios for all performance measures |
| `-l` | `--log` | Print only log of normalizing constant as double |
| `-e` | `--ex` | Print exact normalizing constant numerator and denominator |
| `-g` | `--nc` | Print normalizing constant as double |
| `-t` | `--tput` | Print only throughputs, one per row |
| `-q` | `--qlen` | Print only queue lengths, one per row |
| `-h` | `--help` | Print help message |

#### MOM Solver Options

```bash
./bin/mom [options] model.qn
```

| Option | Long Form | Description |
|--------|-----------|-------------|
| `-v` | `--verbose` | Print exact ratios for all performance measures |
| `-l` | `--log` | Print only log of normalizing constant as double |
| `-e` | `--ex` | Print exact normalizing constant numerator and denominator |
| `-g` | `--nc` | Print normalizing constant as double |
| `-t` | `--tput` | Print only throughputs, one per row |
| `-q` | `--qlen` | Print only queue lengths, one per row |
| `-d` | `--debug` | Enable debug output (progress messages, timing, etc.) |
| `-p digit` | | Apply perturbation at the specified digit (e.g., `-p 5`) |
| `-s seed` | `--seed` | Random seed for perturbation (default: 23000) |
| `-h` | `--help` | Print help message |

#### Examples

```bash
# Get only throughputs
./bin/mva -t models/02_bottleneck_study.qn

# Get only queue lengths  
./bin/mom -q models/02_bottleneck_study.qn

# Get exact normalizing constant (numerator and denominator)
./bin/mva -e models/02_bottleneck_study.qn

# Run MOM solver with debug output
./bin/mom -d models/02_bottleneck_study.qn

# Apply perturbation for approximate solution (MOM only)
./bin/mom -p 5 models/02_bottleneck_study.qn

# Apply perturbation with specific random seed (MOM only)
./bin/mom -p 5 -s 12345 models/02_bottleneck_study.qn
```

**Note**: The perturbation option (`-p`) is only available in the MOM solver and is useful when the exact solution fails due to singular systems. It introduces small randomized numerical perturbations to enable approximate solutions. When perturbation is applied, model parameters are displayed as original integer values with separate perturbation annotations (e.g., `50 eps=1.0e-05`). The `-s` option allows you to specify a random seed for reproducible perturbations.

## References

[1]: Reiser & Lavenberg (1980), *Mean-Value Analysis of Closed Multichain Queuing Networks,* Journal of the ACM 27(2).

[2]: Casale (2006), *An efficient algorithm for the exact analysis of multiclass queueing networks with large population sizes,* Proc. of ACM SIGMETRICS 2006.

## License

This software is provided as-is for academic and research purposes and released as BSD-3.
