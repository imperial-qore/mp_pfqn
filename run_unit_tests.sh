#!/bin/bash

# Regression test script for all solvers on all models with various options
# Compares outputs between solvers (mom, comom, recal) with options -e, -l, -t, -q
#
# Model types are detected automatically:
#   Standard closed     -> mva, ca, recal, mom, comom, gld, mvamx, mvaldmx
#   Closed + MU         -> gld, comomld
#   Mixed (LAMBDA)      -> mvamx
#   Mixed + MU          -> mvaldmx
#
# Usage:
#   ./run_unit_tests.sh                          # Run all models with all options
#   ./run_unit_tests.sh models/03_think.qn       # Run one model with all options
#   ./run_unit_tests.sh models/03_think.qn -- -t # Run one model with -t only

# Set up trap to handle CTRL+C gracefully
trap 'echo -e "\nTest interrupted by user"; exit 130' INT TERM

# cd to the directory containing this script so relative paths work
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd "$SCRIPT_DIR" || exit 1

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
MODELS_DIR="models"
TIMEOUT=120  # Timeout in seconds for each test

# Solver groups
CLOSED_SOLVERS=("./bin/mom" "./bin/mva" "./bin/ca" "./bin/comom" "./bin/recal")
CLOSED_OPTIONS=("-e" "-l" "-t" "-q")

# New solvers that also handle standard closed models (for parity checking)
GLD_SOLVER="./bin/gld"
MVAMX_SOLVER="./bin/mvamx"
MVALDMX_SOLVER="./bin/mvaldmx"
COMOMLD_SOLVER="./bin/comomld"

# Log file
ERROR_LOG="test_errors_$(date +%Y%m%d_%H%M%S).log"

# Function to detect model type by checking for LAMBDA/MU keywords
detect_model_type() {
    local model=$1
    local has_lambda=0
    local has_mu=0
    grep -q "^LAMBDA" "$model" && has_lambda=1
    grep -q "^MU" "$model" && has_mu=1

    if [ $has_lambda -eq 1 ] && [ $has_mu -eq 1 ]; then
        echo "mixed_ld"
    elif [ $has_lambda -eq 1 ]; then
        echo "mixed"
    elif [ $has_mu -eq 1 ]; then
        echo "closed_ld"
    else
        echo "closed"
    fi
}

# Function to detect if a model has multi-server stations (mi > 1)
# GLD and MVA-LD-MX use load-dependent rates instead of station expansion,
# so their normalizing constants and metrics differ on multi-server models.
has_multiserver() {
    local model=$1
    local M=$(sed -n '4p' "$model")
    for i in $(seq 5 $((4 + M))); do
        local mi=$(sed -n "${i}p" "$model" | awk '{print $1}')
        if [ "$mi" -gt 1 ] 2>/dev/null; then
            return 0  # true
        fi
    done
    return 1  # false
}

# Function to get solvers for a model type and option
get_solvers_for_model() {
    local model_type=$1
    local option=$2

    case "$model_type" in
        closed)
            # Standard closed: all original solvers + gld(-e,-l,-g) + mvamx + mvaldmx
            echo "${CLOSED_SOLVERS[*]}"
            ;;
        closed_ld)
            # Closed with MU: gld + comomld (comomld only for -e,-l,-g,-q)
            echo "$GLD_SOLVER $COMOMLD_SOLVER"
            ;;
        mixed)
            # Mixed open/closed: mvamx only
            echo "$MVAMX_SOLVER"
            ;;
        mixed_ld)
            # Mixed + load-dependent: mvaldmx only
            echo "$MVALDMX_SOLVER"
            ;;
    esac
}

# Function to get the subset of options a solver supports
solver_supports_option() {
    local solver=$1
    local option=$2

    case "$(basename $solver)" in
        gld)
            # GLD only supports -e, -l, -g
            [[ "$option" == "-e" || "$option" == "-l" || "$option" == "-g" ]]
            ;;
        comomld)
            # CoMoM-LD supports -e, -l, -g, -q
            [[ "$option" == "-e" || "$option" == "-l" || "$option" == "-g" || "$option" == "-q" ]]
            ;;
        mvamx)
            # MVA-MX supports all standard options
            [[ "$option" == "-e" || "$option" == "-l" || "$option" == "-t" || "$option" == "-q" || "$option" == "-g" || "$option" == "-d" ]]
            ;;
        mvaldmx)
            # MVA-LD-MX supports -t, -q, -d (no normalizing constant)
            [[ "$option" == "-t" || "$option" == "-q" || "$option" == "-d" ]]
            ;;
        *)
            # Original solvers support all standard options
            return 0
            ;;
    esac
}

# Function to print test header
print_header() {
    echo "========================================"
    echo "Regression Tests for MP_PFQN Solvers"
    echo "========================================"
    echo "Date: $(date)"
    echo "Models directory: $MODELS_DIR"
    echo "Timeout: ${TIMEOUT}s"
    echo -e "Error log: $ERROR_LOG"
    echo "========================================"
    echo
}

# Function to run regression tests (compare outputs between solvers)
run_regression_test() {
    local model=$1
    local option=$2
    local model_type=$(detect_model_type "$model")

    echo -e "\n${YELLOW}[${model_type}] $(basename $model) with $option${NC}"

    # Get solvers for this model type
    local solvers_str=$(get_solvers_for_model "$model_type" "$option")
    local solvers=($solvers_str)

    # Store outputs from each solver
    local outputs=()
    local solver_names=()

    for solver in "${solvers[@]}"; do
        local solver_name=$(basename $solver)

        # Check solver exists
        if [ ! -f "$solver" ]; then
            continue
        fi

        # Check solver supports this option
        if ! solver_supports_option "$solver" "$option"; then
            continue
        fi

        timeout $TIMEOUT $solver $option $model > /tmp/regression_${solver_name}_$$.txt 2>/dev/null
        local exit_code=$?
        if [ $exit_code -eq 0 ]; then
            outputs+=("/tmp/regression_${solver_name}_$$.txt")
            solver_names+=("$solver_name")
        elif [ $exit_code -eq 124 ]; then
            echo -e "${RED}TIMEOUT: ${solver_name}${NC}"
            echo "[TIMEOUT] ${solver_name} on $model with $option" >> "$ERROR_LOG"
        else
            echo -e "${RED}ERROR (exit $exit_code): ${solver_name}${NC}"
            echo "[ERROR] ${solver_name} on $model with $option (exit $exit_code)" >> "$ERROR_LOG"
        fi
    done

    # For standard closed models, also run new solvers for cross-validation.
    # GLD and MVA-LD-MX use load-dependent rates (mu) instead of station expansion,
    # so they only match original solvers on single-server models (all mi=1).
    # MVA-MX uses station expansion internally so it always matches.
    if [ "$model_type" == "closed" ]; then
        local is_multiserver=false
        has_multiserver "$model" && is_multiserver=true

        # GLD: cross-validate -e/-l with existing solvers (single-server only)
        if ! $is_multiserver && [ -f "$GLD_SOLVER" ] && solver_supports_option "$GLD_SOLVER" "$option"; then
            timeout $TIMEOUT $GLD_SOLVER $option $model > /tmp/regression_gld_$$.txt 2>/dev/null
            if [ $? -eq 0 ]; then
                outputs+=("/tmp/regression_gld_$$.txt")
                solver_names+=("gld")
            fi
        fi
        # MVA-MX: cross-validate with existing solvers (always matches, uses expansion)
        if [ -f "$MVAMX_SOLVER" ] && solver_supports_option "$MVAMX_SOLVER" "$option"; then
            timeout $TIMEOUT $MVAMX_SOLVER $option $model > /tmp/regression_mvamx_$$.txt 2>/dev/null
            if [ $? -eq 0 ]; then
                outputs+=("/tmp/regression_mvamx_$$.txt")
                solver_names+=("mvamx")
            fi
        fi
        # MVA-LD-MX: cross-validate -t, -q with existing solvers (single-server only)
        if ! $is_multiserver && [ -f "$MVALDMX_SOLVER" ] && solver_supports_option "$MVALDMX_SOLVER" "$option"; then
            timeout $TIMEOUT $MVALDMX_SOLVER $option $model > /tmp/regression_mvaldmx_$$.txt 2>/dev/null
            if [ $? -eq 0 ]; then
                outputs+=("/tmp/regression_mvaldmx_$$.txt")
                solver_names+=("mvaldmx")
            fi
        fi
    fi

    # Compare outputs if we have at least 2
    if [ ${#outputs[@]} -ge 2 ]; then
        local all_match=true
        for ((i=1; i<${#outputs[@]}; i++)); do
            if ! diff -q "${outputs[0]}" "${outputs[$i]}" > /dev/null 2>&1; then
                all_match=false
                echo -e "${RED}Output mismatch:${NC} ${solver_names[0]} vs ${solver_names[$i]}"
                echo "[REGRESSION] Output mismatch for $model with $option: ${solver_names[0]} vs ${solver_names[$i]}" >> "$ERROR_LOG"
                diff "${outputs[0]}" "${outputs[$i]}" >> "$ERROR_LOG" 2>&1
            fi
        done

        if $all_match; then
            echo -e "${GREEN}All ${#outputs[@]} solvers match${NC}"
        fi
    elif [ ${#outputs[@]} -eq 1 ]; then
        echo -e "${GREEN}${solver_names[0]} completed successfully (single solver)${NC}"
    else
        echo -e "${RED}No solvers produced output${NC}"
        echo "[NO OUTPUT] $model with $option" >> "$ERROR_LOG"
    fi

    # Clean up
    rm -f /tmp/regression_*_$$.txt
}

# Parse command-line arguments: [model.qn ...] [-- -opt1 -opt2 ...]
# If models are given, only test those. If options after -- are given, only test those.
parse_args() {
    FILTER_MODELS=()
    FILTER_OPTIONS=()
    local past_separator=false

    for arg in "$@"; do
        if [[ "$arg" == "--" ]]; then
            past_separator=true
            continue
        fi
        if $past_separator; then
            FILTER_OPTIONS+=("$arg")
        else
            FILTER_MODELS+=("$arg")
        fi
    done
}

# Main test execution
main() {
    parse_args "$@"

    # Check if at least some binaries exist
    local found_solver=false
    for solver in "${CLOSED_SOLVERS[@]}" "$GLD_SOLVER" "$COMOMLD_SOLVER" "$MVAMX_SOLVER" "$MVALDMX_SOLVER"; do
        if [ -f "$solver" ]; then
            found_solver=true
            break
        fi
    done
    if ! $found_solver; then
        echo -e "${RED}Error: No solver binaries found. Please run 'make' first.${NC}"
        exit 1
    fi

    # Check if models directory exists
    if [ ! -d "$MODELS_DIR" ]; then
        echo -e "${RED}Error: Models directory $MODELS_DIR not found.${NC}"
        exit 1
    fi

    # Initialize error log
    echo "Regression test errors for test run at $(date)" > "$ERROR_LOG"

    # Build model list: use filter or discover all
    if [ ${#FILTER_MODELS[@]} -gt 0 ]; then
        models=()
        for m in "${FILTER_MODELS[@]}"; do
            # Allow bare filename or path
            if [ -f "$m" ]; then
                models+=("$m")
            elif [ -f "$MODELS_DIR/$m" ]; then
                models+=("$MODELS_DIR/$m")
            else
                echo -e "${RED}Error: Model $m not found.${NC}"
                exit 1
            fi
        done
    else
        models=($(ls $MODELS_DIR/*.qn 2>/dev/null | sort))
    fi

    if [ ${#models[@]} -eq 0 ]; then
        echo -e "${RED}Error: No .qn model files found${NC}"
        exit 1
    fi

    # Build options list: use filter or default
    if [ ${#FILTER_OPTIONS[@]} -gt 0 ]; then
        options=("${FILTER_OPTIONS[@]}")
    else
        options=("${CLOSED_OPTIONS[@]}")
    fi

    # Print header
    print_header
    echo "Found ${#models[@]} models to test"

    # Classify models
    for model in "${models[@]}"; do
        local mtype=$(detect_model_type "$model")
        echo "  $(basename $model) -> $mtype"
    done

    # Run regression tests
    echo -e "\n${YELLOW}Running Regression Tests${NC}"
    echo "========================================"

    for model in "${models[@]}"; do
        for option in "${options[@]}"; do
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
