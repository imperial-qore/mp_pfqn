# Product-Form Queueing Networks (PFQN) Solvers

The mp_pfqn library offers fast C solvers for product-form queueing networks using the GNU Multiple Precision Arithmetic Library (GMP). Exact arithmetic in GMP allows to compute normalizing constants for the equilibrium state probabilities on large closed multiclass models. 

- Current solvers:
  - Mean Value Analysis (MVA) [1]
  - Method of Moments (MoM) [2]
  - Recursion by Chain Algorithm (RECAL) [3]
  - Class-Oriented Method of Moments (CoMoM) [4]

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
./bin/recal models/02_bottleneck_study.qn
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
./bin/recal model.qn    # Recursion by Chain Algorithm
./bin/comom model.qn    # Class-Oriented Method of Moments
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
| `-h` | `--help` | Print help message |


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
| `-d` | `--debug` | Enable debug output (progress messages, timing, etc.) |
| `-p digit` | | Apply perturbation at the specified digit (e.g., `-p 5`) |
| `-s seed` | `--seed` | Random seed for perturbation (default: 23000) |
| `-h` | `--help` | Print help message |

#### CoMoM Solver Options

```bash
./bin/comom [options] models/02_bottleneck_study.qn
```

CoMoM supports the same command-line options as MoM (see table above). The Class-Oriented Method of Moments provides an efficient algorithm specifically designed for multiclass models with large population sizes.

#### Examples

```bash
# Get only throughputs
./bin/mom -t models/02_bottleneck_study.qn

# Get only queue lengths  
./bin/mom -q models/02_bottleneck_study.qn

# Get exact normalizing constant (numerator and denominator)
./bin/mom -e models/02_bottleneck_study.qn

# Run MoM solver with debug output
./bin/mom -d models/02_bottleneck_study.qn

# Apply perturbation for approximate solution (MoM/CoMoM)
./bin/mom -p 5 models/02_bottleneck_study.qn
./bin/comom -p 5 models/02_bottleneck_study.qn

# Apply perturbation with specific random seed (MoM/CoMoM)
./bin/mom -p 5 -s 12345 models/02_bottleneck_study.qn
./bin/comom -p 5 -s 12345 models/02_bottleneck_study.qn
```
**Note**: The perturbation option (`-p`) is available in both the MoM and CoMoM solvers and is useful when the exact solution fails due to singular systems. It introduces small randomized numerical perturbations to enable approximate solutions. When perturbation is applied, model parameters are displayed as original integer values with separate perturbation annotations (e.g., `50 eps=1.0e-05`). The `-s` option allows you to specify a random seed for reproducible perturbations.

## References

[1]: Reiser & Lavenberg (1980), *Mean-Value Analysis of Closed Multichain Queuing Networks,* Journal of the ACM 27(2).

[2]: Casale (2006), *An efficient algorithm for the exact analysis of multiclass queueing networks with large population sizes,* Proc. of ACM SIGMETRICS 2006.

[3]: Conway, A. E. and Georganas, N. D. (1986), *RECAL—a new efficient algorithm for the exact analysis of multiple-chain closed queuing networks,* J. ACM 33, 4 (Oct. 1986), 768–791.

[4]: Casale, G. (2009), *CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models,* IEEE Trans. Softw. Eng. 35, 2 (March 2009), 162–177.

## License

This software is provided as-is for academic and research purposes and released as BSD-3.
