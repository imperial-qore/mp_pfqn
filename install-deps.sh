#!/bin/bash

set -e

echo "=== PFQN Dependencies Installation Script ==="

# Check for sudo privileges upfront
if ! sudo -n true 2>/dev/null; then
    echo "ERROR: This script requires sudo privileges to install system dependencies."
    echo "Please run with: sudo make install-deps"
    echo "Or run: sudo ./install-deps.sh"
    exit 1
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Detect OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
else
    print_error "Unsupported OS: $OSTYPE"
    exit 1
fi

print_status "Detected OS: $OS"

# Install system dependencies
print_status "Installing system dependencies..."

if [[ "$OS" == "linux" ]]; then
    if command -v apt-get &> /dev/null; then
        sudo apt-get update
        sudo apt-get install -y build-essential autotools-dev autoconf automake libtool m4 cmake pkg-config git wget \
                               libblas-dev liblapack-dev libopenblas-dev gfortran
    elif command -v yum &> /dev/null; then
        sudo yum install -y gcc gcc-c++ gcc-gfortran make autoconf automake libtool m4 cmake pkgconfig git wget \
                           blas-devel lapack-devel openblas-devel
    elif command -v dnf &> /dev/null; then
        sudo dnf install -y gcc gcc-c++ gcc-gfortran make autoconf automake libtool m4 cmake pkgconfig git wget \
                           blas-devel lapack-devel openblas-devel
    else
        print_error "Unsupported Linux distribution. Please install build tools manually."
        exit 1
    fi
elif [[ "$OS" == "macos" ]]; then
    if ! command -v brew &> /dev/null; then
        print_error "Homebrew not found. Please install Homebrew first: https://brew.sh"
        exit 1
    fi
    brew install autoconf automake libtool cmake pkg-config git wget openblas lapack
fi

# Create deps directory
DEPS_DIR="$(pwd)/deps"
mkdir -p "$DEPS_DIR"
cd "$DEPS_DIR"

print_status "Dependencies will be installed in: $DEPS_DIR"

print_status "Installing GMP..."

# Install GMP
if [[ ! -f "gmp-6.3.0.tar.xz" ]]; then
    wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
fi

if [[ ! -d "gmp-6.3.0" ]]; then
    tar -xf gmp-6.3.0.tar.xz
fi

cd gmp-6.3.0
if [[ ! -f "$DEPS_DIR/local/lib/libgmp.so" ]] && [[ ! -f "$DEPS_DIR/local/lib/libgmp.a" ]]; then
    # Clean any previous partial builds
    make distclean 2>/dev/null || true
    
    # Configure with absolute path resolution
    PREFIX_PATH=$(cd "$DEPS_DIR/local" && pwd) || PREFIX_PATH="$DEPS_DIR/local"
    ./configure --prefix="$PREFIX_PATH" --enable-cxx --disable-shared
    make -j$(nproc 2>/dev/null || echo 4)
    make install
fi
cd "$DEPS_DIR"

print_status "Installing MPFR..."

# Install MPFR
if [[ ! -f "mpfr-4.2.1.tar.xz" ]]; then
    wget https://www.mpfr.org/mpfr-4.2.1/mpfr-4.2.1.tar.xz || \
    wget https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz || \
    wget https://ftpmirror.gnu.org/mpfr/mpfr-4.2.1.tar.xz
fi

if [[ ! -d "mpfr-4.2.1" ]]; then
    tar -xf mpfr-4.2.1.tar.xz
fi

cd mpfr-4.2.1
if [[ ! -f "$DEPS_DIR/local/lib/libmpfr.so" ]] && [[ ! -f "$DEPS_DIR/local/lib/libmpfr.a" ]]; then
    # Clean any previous partial builds
    make distclean 2>/dev/null || true
    
    # Configure with absolute path resolution
    PREFIX_PATH=$(cd "$DEPS_DIR/local" && pwd) || PREFIX_PATH="$DEPS_DIR/local"
    ./configure --prefix="$PREFIX_PATH" --with-gmp="$PREFIX_PATH" --disable-shared
    make -j$(nproc 2>/dev/null || echo 4)
    make install
fi
cd "$DEPS_DIR"

# Install Givaro (required by fflas-ffpack)
print_status "Installing Givaro..."

if [[ ! -d "givaro" ]]; then
    git clone https://github.com/linbox-team/givaro.git
fi

cd givaro
PREFIX_PATH=$(cd "$DEPS_DIR/local" && pwd) || PREFIX_PATH="$DEPS_DIR/local"
export PKG_CONFIG_PATH="$PREFIX_PATH/lib/pkgconfig:$PKG_CONFIG_PATH"
export CPPFLAGS="-I$PREFIX_PATH/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX_PATH/lib $LDFLAGS"

if [[ ! -f "$PREFIX_PATH/lib/libgivaro.a" ]]; then
    if [[ ! -f "configure" ]]; then
        ./autogen.sh
    fi
    
    ./configure --prefix="$PREFIX_PATH" \
               --with-gmp="$PREFIX_PATH" \
               --disable-shared
    
    make -j$(nproc 2>/dev/null || echo 4)
    make install
fi

cd "$DEPS_DIR"

# Install fflas-ffpack (required by LinBox)
print_status "Installing fflas-ffpack..."

if [[ ! -d "fflas-ffpack" ]]; then
    git clone https://github.com/linbox-team/fflas-ffpack.git
fi

cd fflas-ffpack

if [[ ! -f "$PREFIX_PATH/lib/pkgconfig/fflas-ffpack.pc" ]]; then
    if [[ ! -f "configure" ]]; then
        ./autogen.sh
    fi
    
    ./configure --prefix="$PREFIX_PATH" \
               --with-gmp="$PREFIX_PATH" \
               --with-givaro="$PREFIX_PATH" \
               --with-blas-libs="-lopenblas" \
               --disable-shared
    
    make -j$(nproc 2>/dev/null || echo 4)
    make install
fi

cd "$DEPS_DIR"

# Install LinBox
print_status "Installing LinBox..."

if [[ ! -d "linbox" ]]; then
    git clone https://github.com/linbox-team/linbox.git
fi

cd linbox

if [[ ! -f "$PREFIX_PATH/lib/liblinbox.a" ]]; then
    if [[ ! -f "configure" ]]; then
        ./autogen.sh
    fi
    
    ./configure --prefix="$PREFIX_PATH" \
               --with-gmp="$PREFIX_PATH" \
               --with-mpfr="$PREFIX_PATH" \
               --with-givaro="$PREFIX_PATH" \
               --with-fflas-ffpack="$PREFIX_PATH" \
               --with-blas-libs="-lopenblas" \
               --disable-shared
    
    make -j$(nproc 2>/dev/null || echo 4)
    make install
fi

cd "$DEPS_DIR/.."

# Create environment setup script
print_status "Creating environment setup script..."

cat > setup-env.sh << 'EOF'
#!/bin/bash
# Source this script to set up the environment for building PFQN

DEPS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/deps" && pwd)"

if [[ -d "$DEPS_DIR/local" ]]; then
    export PKG_CONFIG_PATH="$DEPS_DIR/local/lib/pkgconfig:$PKG_CONFIG_PATH"
    export CPPFLAGS="-I$DEPS_DIR/local/include $CPPFLAGS"
    export LDFLAGS="-L$DEPS_DIR/local/lib $LDFLAGS"
    export LD_LIBRARY_PATH="$DEPS_DIR/local/lib:$LD_LIBRARY_PATH"
    echo "Using dependencies from $DEPS_DIR/local"
else
    echo "No local dependencies found, using system libraries"
fi

echo "Environment configured for PFQN build"
EOF

chmod +x setup-env.sh

print_status "Installation complete!"
print_status "Dependencies installed in: $DEPS_DIR/local"
print_status "To build the project:"
print_status "  make"
print_status ""
print_status "Or if you need to set up environment manually:"
print_status "  source ./setup-env.sh"
print_status "  make"