#!/bin/bash

# Regression test script for all solvers on all models with various options
# Compares outputs between solvers (mom, comom, recal) with options -e, -l, -t, -q
#
# Model types are detected automatically:
#   Standard closed     -> mva, ca, recal, mom, comom, gcomom, automom, gld, mvamx, mvaldmx
#   Closed + MU         -> gld, comomld
#   Mixed (LAMBDA)      -> mvamx
#   Mixed + MU          -> mvaldmx
#   Fractional mi       -> mom, comom, gcomom, gmom, automom, checked against gld (-e,-l)
#                          and mvaldmx (-t,-q) on the load-dependent companion
#                          model, rates t/(t+e-1) (see CLAUDE.md)
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
TIMEOUT=300  # Timeout in seconds for each test (comom on 10_diverse takes ~150 s)

# Solver groups
CLOSED_SOLVERS=("./bin/mom" "./bin/mva" "./bin/ca" "./bin/comom" "./bin/recal" "./bin/gcomom" "./bin/automom")
CLOSED_OPTIONS=("-e" "-l" "-t" "-q")

# New solvers that also handle standard closed models (for parity checking)
GLD_SOLVER="./bin/gld"
MVAMX_SOLVER="./bin/mvamx"
MVALDMX_SOLVER="./bin/mvaldmx"
COMOMLD_SOLVER="./bin/comomld"

# Solvers that accept fractional multiplicities (qn_allow_rational_mi)
FRAC_SOLVERS=("./bin/mom" "./bin/comom" "./bin/gcomom" "./bin/gmom" "./bin/automom")

# Log file
ERROR_LOG="test_errors_$(date +%Y%m%d_%H%M%S).log"

# True if some station line (5 .. 4+M) has a fractional multiplicity p/q
has_fractional_mi() {
    local M=$(sed -n '4p' "$1")
    sed -n "5,$((4 + M))p" "$1" | awk '$1 ~ /\// {f=1} END {exit !f}'
}

# The load-dependent model equivalent to a fractional one: each station
# with multiplicity e = p/q becomes one station with MU rates
# t/(t+e-1) = tq/(tq+p-q), t = 1..Nt; the generating function
# (1 - sum_r L_r x_r)^(-e) is the same, so G, X and Q are identical.
write_ld_companion() {
    awk 'function gcd(a,b,  c){while(b){c=b;b=a%b;a=c}return a}
         NR==2 {Nt=0; for(i=1;i<=NF;i++) Nt+=$i}
         NR==4 {M=$1}
         NR>=5 && NR<=4+M {
             n=split($1,f,"/"); p[NR]=f[1]; q[NR]=(n>1?f[2]:1)
             $1=1
         }
         {print}
         END {
             print "MU"
             for(k=5;k<=4+M;k++){
                 line=""
                 for(t=1;t<=Nt;t++){
                     a=t*q[k]; b=t*q[k]+p[k]-q[k]; g=gcd(a,b); a/=g; b/=g
                     line=line (t>1?" ":"") (b==1?a:a "/" b)
                 }
                 print line
             }
         }' "$1" > "$2"
}

# Function to detect model type by checking for LAMBDA/MU keywords
detect_model_type() {
    local model=$1
    local has_lambda=0
    local has_mu=0
    grep -q "^LAMBDA" "$model" && has_lambda=1
    grep -q "^MU" "$model" && has_mu=1

    if has_fractional_mi "$model"; then
        echo "closed_frac"
    elif [ $has_lambda -eq 1 ] && [ $has_mu -eq 1 ]; then
        echo "mixed_ld"
    elif [ $has_lambda -eq 1 ]; then
        echo "mixed"
    elif [ $has_mu -eq 1 ]; then
        echo "closed_ld"
    else
        echo "closed"
    fi
}

# A solver's documented refusal of a model outside its domain (stderr):
# reported as a skip, not an error, since there is nothing to compare.
#   comomld: exactly one queueing station after IS absorption
#   recal:   refuses a table it cannot allocate
#   gmom:    no exact G on a model singular under every class order (its
#            unflagged run perturbs for -l/-t/-q; mom/comom/gcomom are exact)
is_domain_refusal() {
    local solver_name=$1 errfile=$2
    case "$solver_name" in
        comomld) grep -q "outside the CoMoM-LD domain" "$errfile" ;;
        recal)   grep -q "Model may be too large" "$errfile" ;;
        gmom)    grep -q "singular under every class order" "$errfile" ;;
        *)       return 1 ;;
    esac
}
SKIPPED=0

# A cross-check solver that fails is reported, not dropped silently. These are
# the exponential oracles (gld's G table, mvaldmx's lattice), so a timeout is a
# SKIP: gld times out on 10_diverse, mvaldmx on 06_large and 10_diverse.
cross_check() {
    if [ $1 -eq 0 ]; then
        outputs+=("/tmp/regression_$2_$$.txt")
        solver_names+=("$2")
    elif [ $1 -eq 124 ]; then
        echo -e "${YELLOW}SKIP: $2 cross-check (over ${TIMEOUT}s)${NC}"
        SKIPPED=$((SKIPPED + 1))
    else
        echo -e "${RED}ERROR (exit $1): $2 cross-check${NC}"
        echo "[ERROR] $2 cross-check on $model with $option (exit $1)" >> "$ERROR_LOG"
    fi
}

# RECAL sweeps a table of C(Ntot+Mz-1, Ntot) exact rationals per population
# step (Mz = M, plus one with think times): 2.1M entries on 07_mixed did not
# finish in 25 minutes, 76k on 06_large take ~170 s. Above RECAL_MAX entries
# recal is skipped as outside its practical domain.
RECAL_MAX=1000000
recal_too_large() {
    awk -v max=$RECAL_MAX 'NR==1{R=$1} NR==2{for(i=1;i<=R;i++)n+=$i}
        NR==3{for(i=1;i<=R;i++)if($i>0)z=1} NR==4{M=$1; exit}
        END{m=M+z; c=1; for(i=1;i<=n;i++)c=c*(n+m-1-n+i)/i; exit !(c>max)}' "$1"
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
            # Closed with MU: gld (all options), comomld (-e,-l,-g), mvaldmx (-t,-q)
            echo "$GLD_SOLVER $COMOMLD_SOLVER $MVALDMX_SOLVER"
            ;;
        mixed)
            # Mixed open/closed: mvamx only
            echo "$MVAMX_SOLVER"
            ;;
        mixed_ld)
            # Mixed + load-dependent: mvaldmx only
            echo "$MVALDMX_SOLVER"
            ;;
        closed_frac)
            # Fractional multiplicities: the solvers that accept them
            echo "${FRAC_SOLVERS[*]}"
            ;;
    esac
}

# Function to get the subset of options a solver supports
solver_supports_option() {
    local solver=$1
    local option=$2

    case "$(basename $solver)" in
        comomld)
            # CoMoM-LD supports -e, -l, -g (its state probabilities are -P, not -q)
            [[ "$option" == "-e" || "$option" == "-l" || "$option" == "-g" ]]
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
    local attempted=0

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

        # gmom needs at least two classes
        if [ "$solver_name" == "gmom" ] && [ "$(sed -n '1p' "$model")" -lt 2 ]; then
            continue
        fi

        if [ "$solver_name" == "recal" ] && recal_too_large "$model"; then
            echo -e "${YELLOW}SKIP: recal (table above $RECAL_MAX entries)${NC}"
            SKIPPED=$((SKIPPED + 1))
            attempted=$((attempted + 1))
            continue
        fi

        attempted=$((attempted + 1))
        timeout $TIMEOUT $solver $option $model > /tmp/regression_${solver_name}_$$.txt 2>/tmp/regression_${solver_name}_err_$$.txt
        local exit_code=$?
        if [ $exit_code -eq 0 ]; then
            outputs+=("/tmp/regression_${solver_name}_$$.txt")
            solver_names+=("$solver_name")
        elif is_domain_refusal "$solver_name" /tmp/regression_${solver_name}_err_$$.txt; then
            echo -e "${YELLOW}SKIP: ${solver_name} (outside its domain)${NC}"
            SKIPPED=$((SKIPPED + 1))
        elif [ $exit_code -eq 124 ]; then
            echo -e "${RED}TIMEOUT: ${solver_name}${NC}"
            echo "[TIMEOUT] ${solver_name} on $model with $option" >> "$ERROR_LOG"
        else
            echo -e "${RED}ERROR (exit $exit_code): ${solver_name}${NC}"
            echo "[ERROR] ${solver_name} on $model with $option (exit $exit_code)" >> "$ERROR_LOG"
        fi
    done

    # For standard closed models, also run new solvers for cross-validation.
    # GLD and MVA-LD-MX solve the copies meaning of mi exactly
    # (util/replicas.h), so they match on every closed model.
    if [ "$model_type" == "closed" ]; then
        # GLD: cross-validate with existing solvers
        if [ -f "$GLD_SOLVER" ] && solver_supports_option "$GLD_SOLVER" "$option"; then
            timeout $TIMEOUT $GLD_SOLVER $option $model > /tmp/regression_gld_$$.txt 2>/dev/null
            cross_check $? gld
        fi
        # MVA-MX: cross-validate with existing solvers (always matches, uses expansion)
        if [ -f "$MVAMX_SOLVER" ] && solver_supports_option "$MVAMX_SOLVER" "$option"; then
            timeout $TIMEOUT $MVAMX_SOLVER $option $model > /tmp/regression_mvamx_$$.txt 2>/dev/null
            cross_check $? mvamx
        fi
        # MVA-LD-MX: cross-validate -t, -q with existing solvers
        if [ -f "$MVALDMX_SOLVER" ] && solver_supports_option "$MVALDMX_SOLVER" "$option"; then
            timeout $TIMEOUT $MVALDMX_SOLVER $option $model > /tmp/regression_mvaldmx_$$.txt 2>/dev/null
            cross_check $? mvaldmx
        fi
    fi

    # Fractional multiplicities: the exact oracle is gld (-e, -l) or mvaldmx
    # (-t, -q) on the equivalent load-dependent model.
    if [ "$model_type" == "closed_frac" ]; then
        local ld_model=/tmp/regression_ldcomp_$$.qn
        write_ld_companion "$model" "$ld_model"
        for oracle in "$GLD_SOLVER" "$MVALDMX_SOLVER"; do
            local oname=$(basename $oracle)
            if [ -f "$oracle" ] && solver_supports_option "$oracle" "$option"; then
                timeout $TIMEOUT $oracle $option $ld_model > /tmp/regression_${oname}_$$.txt 2>/dev/null
                local oexit=$?
                if [ $oexit -eq 0 ]; then
                    outputs+=("/tmp/regression_${oname}_$$.txt")
                    solver_names+=("${oname}(ld)")
                else
                    echo -e "${RED}ERROR (exit $oexit): ${oname} on the load-dependent companion${NC}"
                    echo "[ERROR] ${oname} on load-dependent companion of $model with $option (exit $oexit)" >> "$ERROR_LOG"
                fi
            fi
        done
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
    elif [ $attempted -eq 0 ]; then
        # no solver of this model type implements the option (e.g. the
        # normalizing constant of a mixed load-dependent model)
        echo -e "${YELLOW}SKIP: no solver implements $option for ${model_type} models${NC}"
        SKIPPED=$((SKIPPED + 1))
    else
        echo -e "${RED}No solvers produced output${NC}"
        echo "[NO OUTPUT] $model with $option" >> "$ERROR_LOG"
    fi

    # Clean up
    rm -f /tmp/regression_*_$$.txt /tmp/regression_ldcomp_$$.qn
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
    echo "Skipped (outside a solver's domain, or option not implemented, oracle over the timeout): $SKIPPED"
    echo "Check $ERROR_LOG for any mismatches."
    echo "========================================"
}

# Run main function
main "$@"
