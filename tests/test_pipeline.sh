#!/bin/bash
# Test script for the complete RNA scoring pipeline
# Tests all CLI commands: access, extract, train, score, and plot

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test results
declare -A RESULTS

echo ""
echo "========================================================================"
echo "RNA Scoring Pipeline - Full Test Suite"
echo "========================================================================"
echo ""

# Create temporary directory for test outputs
TMPDIR=$(mktemp -d -p /tmp rna_test_XXXXXX)
trap "rm -rf $TMPDIR" EXIT

echo "Using temporary directory: $TMPDIR"
echo ""

# Function to run command and log results
run_test() {
    local test_name="$1"
    local description="$2"
    shift 2
    local cmd=("$@")
    
    echo ""
    echo "========================================================================"
    echo "🔧 $description"
    echo "========================================================================"
    echo "$ ${cmd[*]}"
    echo ""
    
    if "${cmd[@]}"; then
        echo -e "${GREEN} $description - SUCCESS${NC}"
        RESULTS["$test_name"]=1
        return 0
    else
        echo -e "${RED} $description - FAILED${NC}"
        RESULTS["$test_name"]=0
        return 1
    fi
}

# Test 1: List available structures
echo ""
echo "========================================================================"
echo "TEST 1: List available RNA structures (no download)"
echo "========================================================================"
if python3 -m rna_score.cli access --list-only -n 5; then
    echo -e "${GREEN} List structures - SUCCESS${NC}"
    RESULTS["list"]=1
else
    echo -e "${RED} List structures - FAILED${NC}"
    RESULTS["list"]=0
fi

# Test 2: Download structures
DOWNLOAD_DIR="$TMPDIR/structures"
echo ""
echo "========================================================================"
echo "TEST 2: Download RNA structures"
echo "========================================================================"
echo "$ python3 -m rna_score.cli access -n 50 -f cif -o $DOWNLOAD_DIR -w 2"
echo ""

if python3 -m rna_score.cli access -n 100 -f cif -o "$DOWNLOAD_DIR" -w 2; then
    echo -e "${GREEN} Download structures - SUCCESS${NC}"
    RESULTS["access"]=1
    
    # Check what was downloaded
    if [ -d "$DOWNLOAD_DIR/mmcif" ]; then
        CIF_COUNT=$(find "$DOWNLOAD_DIR/mmcif" -name "*.cif" 2>/dev/null | wc -l)
        echo "   Found $CIF_COUNT .cif files"
    fi
else
    echo -e "${RED} Download structures - FAILED${NC}"
    RESULTS["access"]=0
    echo " Access step failed. Skipping downstream tests."
    
    # Print summary and exit
    echo ""
    echo "========================================================================"
    echo "TEST SUMMARY"
    echo "========================================================================"
    for test in list access extract train score plot; do
        if [ "${RESULTS[$test]}" = "1" ]; then
            echo -e "${GREEN} PASS${NC}   $test"
        elif [ -n "${RESULTS[$test]}" ]; then
            echo -e "${RED} FAIL${NC}   $test"
        else
            echo -e "${YELLOW}⊘ SKIP${NC}   $test"
        fi
    done
    exit 1
fi

# Test 3: Extract distances
DIST_DIR="$TMPDIR/distances"
PDB_IDS_FILE="$TMPDIR/pdb_ids.txt"

echo ""
echo "========================================================================"
echo "TEST 3: Extract distances from structures"
echo "========================================================================"

# Prefer using the downloaded mmcif folder to maximize coverage
if [ -d "$DOWNLOAD_DIR/mmcif" ]; then
    echo "$ python3 -m rna_score.cli extract --folder $DOWNLOAD_DIR/mmcif --format mmcif --out-dir $DIST_DIR"
    echo ""

    if python3 -m rna_score.cli extract --folder "$DOWNLOAD_DIR/mmcif" --format mmcif --out-dir "$DIST_DIR"; then
        echo -e "${GREEN} Extract distances - SUCCESS${NC}"
        RESULTS["extract"]=1
        echo "   Histogram files present:"
        ls "$DIST_DIR"/hist_*.csv 2>/dev/null || true
    else
        echo -e "${RED} Extract distances - FAILED${NC}"
        RESULTS["extract"]=0
        echo " Extract step failed."
    fi
else
    # Fallback: use a small list of IDs if no downloads present
    cat > "$PDB_IDS_FILE" << 'EOF'
1J8G
1P79
1Q9A
EOF
    echo "   Using fallback PDB IDs: 1J8G, 1P79, 1Q9A"
    echo "$ python3 -m rna_score.cli extract --list $PDB_IDS_FILE --format mmcif --out-dir $DIST_DIR"
    echo ""

    if python3 -m rna_score.cli extract --list "$PDB_IDS_FILE" --format mmcif --out-dir "$DIST_DIR"; then
        echo -e "${GREEN} Extract distances - SUCCESS${NC}"
        RESULTS["extract"]=1
    else
        echo -e "${RED} Extract distances - FAILED${NC}"
        RESULTS["extract"]=0
        echo " Extract step failed."
    fi
fi

# Early exit on extraction failure
if [ "${RESULTS["extract"]}" != "1" ]; then
    echo ""
    echo "========================================================================"
    echo "TEST SUMMARY"
    echo "========================================================================"
    for test in list access extract train score plot; do
        if [ "${RESULTS[$test]}" = "1" ]; then
            echo -e "${GREEN} PASS${NC}   $test"
        elif [ -n "${RESULTS[$test]}" ]; then
            echo -e "${RED} FAIL${NC}   $test"
        else
            echo -e "${YELLOW}⊘ SKIP${NC}   $test"
        fi
    done
    exit 1
fi

# Test 4: Train scoring function
TRAIN_DIR="$TMPDIR/training"
echo ""
echo "========================================================================"
echo "TEST 4: Train scoring function"
echo "========================================================================"
echo "$ python3 -m rna_score.cli train --input-dir $DIST_DIR --output-dir $TRAIN_DIR"
echo ""

if python3 -m rna_score.cli train --input-dir "$DIST_DIR" --output-dir "$TRAIN_DIR"; then
    echo -e "${GREEN} Train scoring function - SUCCESS${NC}"
    RESULTS["train"]=1
else
    echo -e "${RED} Train scoring function - FAILED${NC}"
    RESULTS["train"]=0
    echo " Train step failed."
    
    # Print partial summary
    echo ""
    echo "========================================================================"
    echo "TEST SUMMARY"
    echo "========================================================================"
    for test in list access extract train score plot; do
        if [ "${RESULTS[$test]}" = "1" ]; then
            echo -e "${GREEN} PASS${NC}   $test"
        elif [ -n "${RESULTS[$test]}" ]; then
            echo -e "${RED} FAIL${NC}   $test"
        else
            echo -e "${YELLOW}⊘ SKIP${NC}   $test"
        fi
    done
    exit 1
fi

# Test 5: Score structures
SCORE_DIR="$TMPDIR/scores"

echo ""
echo "========================================================================"
echo "TEST 5: Score structures"
echo "========================================================================"

# Use downloaded CIF files directory for scoring
if [ -d "$DOWNLOAD_DIR/mmcif" ]; then
    CIF_COUNT=$(find "$DOWNLOAD_DIR/mmcif" -name "*.cif" 2>/dev/null | wc -l)
    echo "   Found $CIF_COUNT CIF files in $DOWNLOAD_DIR/mmcif"
    
    echo "$ python3 -m rna_score.cli score --folder $DOWNLOAD_DIR/mmcif --tables $TRAIN_DIR --format mmcif --output $SCORE_DIR/scores.csv"
    echo ""

    if python3 -m rna_score.cli score --folder "$DOWNLOAD_DIR/mmcif" --tables "$TRAIN_DIR" --format mmcif --output "$SCORE_DIR/scores.csv"; then
        echo -e "${GREEN} Score structures - SUCCESS${NC}"
        RESULTS["score"]=1
    else
        echo -e "${RED} Score structures - FAILED${NC}"
        RESULTS["score"]=0
        echo " Score step failed."
    fi
else
    echo -e "${RED} Score structures - SKIPPED (no CIF directory)${NC}"
    RESULTS["score"]=0
fi

# Test 6: Plot scores
PLOT_DIR="$TMPDIR/plots"
echo ""
echo "========================================================================"
echo "TEST 6: Plot score distributions"
echo "========================================================================"
echo "$ python3 -m rna_score.cli plot --input-dir $TRAIN_DIR --output-dir $PLOT_DIR"
echo ""

if python3 -m rna_score.cli plot --input-dir "$TRAIN_DIR" --output-dir "$PLOT_DIR"; then
    echo -e "${GREEN} Plot scores - SUCCESS${NC}"
    RESULTS["plot"]=1
else
    echo -e "${RED} Plot scores - FAILED${NC}"
    RESULTS["plot"]=0
fi

# Print final summary
echo ""
echo ""
echo "========================================================================"
echo "TEST SUMMARY"
echo "========================================================================"

PASSED=0
FAILED=0
SKIPPED=0

for test in list access extract train score plot; do
    if [ "${RESULTS[$test]}" = "1" ]; then
        echo -e "${GREEN} PASS${NC}   $test"
        ((PASSED++))
    elif [ -n "${RESULTS[$test]}" ]; then
        echo -e "${RED} FAIL${NC}   $test"
        ((FAILED++))
    else
        echo -e "${YELLOW}⊘ SKIP${NC}   $test"
        ((SKIPPED++))
    fi
done

TOTAL=$((PASSED + FAILED))
echo ""
echo "$PASSED/$TOTAL tests passed"

if [ $FAILED -eq 0 ] && [ $PASSED -gt 0 ]; then
    echo ""
    echo -e "${GREEN} All tests passed!${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}  $FAILED test(s) failed${NC}"
    exit 1
fi
