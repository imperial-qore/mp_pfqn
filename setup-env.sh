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
