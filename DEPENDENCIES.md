# Dependencies

This document details all dependencies required for building and running the PFQN solvers.

## Core Dependencies

### Required Libraries

- **GMP 6.3.0**: GNU Multiple Precision Arithmetic Library
  - Used for arbitrary precision arithmetic in all solvers
  - Required for numerical stability in large models
  
- **MPFR 4.2.1**: Multiple Precision Floating-point Library
  - Provides arbitrary precision floating-point operations
  - Built on top of GMP

### Build-time Dependencies

- Standard math library (libm)

### System Dependencies

- **Build Tools**:
  - gcc/g++ (C/C++ compiler)
  - GNU Make
  - autotools (autoconf, automake, libtool)
  - m4 macro processor
  - cmake (for building dependencies)
  - pkg-config
  - git
  - wget

- **Linear Algebra Libraries**:
  - BLAS (Basic Linear Algebra Subprograms)
  - LAPACK (Linear Algebra Package)
  - OpenBLAS (optimized BLAS implementation)
  - gfortran (Fortran compiler for BLAS/LAPACK)

## Installation Options

### Option 1: Automatic Installation (Recommended)

```bash
# Install all dependencies automatically (requires sudo)
sudo make install-deps

# Build the project
make
```

This will automatically install all dependencies, building them from source as static libraries in `deps/local/`.

### Option 2: System Package Installation

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install build-essential autotools-dev autoconf automake libtool m4 \
                     cmake pkg-config git wget gfortran \
                     libgmp-dev libmpfr-dev libblas-dev liblapack-dev libopenblas-dev
```

#### CentOS/RHEL/Fedora
```bash
sudo yum install gcc gcc-c++ gcc-gfortran make autoconf automake libtool m4 \
                 cmake pkgconfig git wget \
                 gmp-devel mpfr-devel blas-devel lapack-devel openblas-devel
```

#### macOS
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install autoconf automake libtool cmake pkg-config git wget \
             gmp mpfr openblas lapack
```


### Option 3: Custom Installation

If you have dependencies installed in non-standard locations:
```bash
# Set environment variables
export CPPFLAGS="-I/path/to/include $CPPFLAGS"
export LDFLAGS="-L/path/to/lib $LDFLAGS"
export PKG_CONFIG_PATH="/path/to/lib/pkgconfig:$PKG_CONFIG_PATH"

# Build
make
```

## Dependency Detection

The build system automatically detects dependencies in this order:
1. **Local installation** in `deps/local/` (from `sudo make install-deps`)
2. **System libraries** (requires manual installation)

You can verify which dependencies are being used:
```bash
make check-deps
```

## Solver Requirements

- **All solvers**: Require GMP and MPFR

## Troubleshooting

### Common Issues

1. **"cannot find -lmpfr" error**: Dependencies not installed. Run `sudo make install-deps`

2. **"BLAS not found" error**: BLAS libraries missing. Install with:
   - Ubuntu/Debian: `sudo apt-get install libopenblas-dev`
   - CentOS/RHEL: `sudo yum install openblas-devel`  
   - macOS: `brew install openblas`

3. **Permission denied**: The install script requires sudo for system packages:
   ```bash
   sudo make install-deps  # Correct
   make install-deps       # Will fail
   ```

4. **WSL path issues**: The install script handles WSL mount points automatically

### Environment Setup

After running `sudo make install-deps`, a `setup-env.sh` script is created. If needed:
```bash
source ./setup-env.sh  # Sets up library paths
make
```

## Version Information

The automatic installation (Option 1) installs specific versions that have been tested:
- GMP 6.3.0
- MPFR 4.2.1

System package versions may vary depending on your distribution.