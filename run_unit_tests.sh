#!/bin/bash

# Regression test script for all solvers on all models with various options
# Compares outputs between solvers (mom, comom, recal) with options -e, -l, -t, -q

# Set up trap to handle CTRL+C gracefully
trap 'echo -e "\nTest interrupted by user"; exit 130' INT TERM

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
MODELS_DIR="models"
SOLVERS=("./bin/mom" "./bin/comom" "./bin/recal")
OPTIONS=("-e" "-l" "-t" "-q")
TIMEOUT=10  # Timeout in seconds for each test

# Log file
ERROR_LOG="test_errors_$(date +%Y%m%d_%H%M%S).log"

# Function to print test header
print_header() {
    echo "========================================"
    echo "Regression Tests for MP_PFQN Solvers"
    echo "========================================"
    echo "Date: $(date)"
    echo "Models directory: $MODELS_DIR"
    echo "Solvers: ${SOLVERS[*]}"
    echo "Options to test: ${OPTIONS[*]}"
    echo "Timeout: ${TIMEOUT}s"
    echo -e "Error log: $ERROR_LOG"
    echo "========================================"
    echo
}

# Function to run regression tests (compare outputs between solvers)
run_regression_test() {
    local model=$1
    local option=$2
    
    echo -e "\n${YELLOW}Regression test for $(basename $model) with option $option${NC}"
    
    # Store outputs from each solver
    local outputs=()
    local solver_names=()
    
    for solver in "${SOLVERS[@]}"; do
        local solver_name=$(basename $solver)
        
        # Skip unsupported options
        if [[ "$solver_name" == "recal" && ("$option" == "-p" || "$option" == "-s") ]]; then
            continue
        fi
        
        timeout $TIMEOUT $solver $option $model > /tmp/regression_${solver_name}_$$.txt 2>/dev/null
        if [ $? -eq 0 ]; then
            outputs+=("/tmp/regression_${solver_name}_$$.txt")
            solver_names+=("$solver_name")
        fi
    done
    
    # Compare outputs if we have at least 2
    if [ ${#outputs[@]} -ge 2 ]; then
        local all_match=true
        for ((i=1; i<${#outputs[@]}; i++)); do
            if ! diff -q "${outputs[0]}" "${outputs[$i]}" > /dev/null 2>&1; then
                all_match=false
                echo -e "${RED}Output mismatch:${NC} ${solver_names[0]} vs ${solver_names[$i]}"
                echo "[REGRESSION] Output mismatch for $model with $option: ${solver_names[0]} vs ${solver_names[$i]}" >> "$ERROR_LOG"
            fi
        done
        
        if $all_match; then
            echo -e "${GREEN}All solvers produce consistent output${NC}"
        fi
    fi
    
    # Clean up
    rm -f /tmp/regression_*_$$.txt
}

# Main test execution
main() {
    # Check if binaries exist
    for solver in "${SOLVERS[@]}"; do
        if [ ! -f "$solver" ]; then
            echo -e "${RED}Error: Solver $solver not found. Please run 'make' first.${NC}"
            exit 1
        fi
    done
    
    # Check if models directory exists
    if [ ! -d "$MODELS_DIR" ]; then
        echo -e "${RED}Error: Models directory $MODELS_DIR not found.${NC}"
        exit 1
    fi
    
    # Initialize error log
    echo "Regression test errors for test run at $(date)" > "$ERROR_LOG"
    
    # Print header
    print_header
    
    # Get all model files
    models=($(ls $MODELS_DIR/*.qn 2>/dev/null | sort))
    
    if [ ${#models[@]} -eq 0 ]; then
        echo -e "${RED}Error: No .qn model files found in $MODELS_DIR${NC}"
        exit 1
    fi
    
    echo "Found ${#models[@]} models to test"
    
    # Run regression tests for key options
    echo -e "\n${YELLOW}Running Regression Tests${NC}"
    echo "========================================"
    
    for model in "${models[@]}"; do
        for option in "${OPTIONS[@]}"; do
            run_regression_test "$model" "$option"
        done
    done
    
    echo -e "\n========================================"
    echo "Regression tests completed."
    echo "Check $ERROR_LOG for any mismatches."
    echo "========================================"
}

# Run main function
main "$@"