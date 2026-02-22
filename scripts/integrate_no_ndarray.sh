#!/bin/bash
# integrate_no_ndarray.sh - Incremental integration of no_ndarray branch
#
# This script helps safely integrate the no_ndarray branch into master
# with careful testing and validation at each step.
#
# Usage: ./scripts/integrate_no_ndarray.sh [step]
#   step 1: Create integration branch and backup
#   step 2: Cherry-pick infrastructure improvements
#   step 3: Integrate GPU-specific code
#   step 4: Update build system
#   step 5: Testing and validation

set -e

INTEGRATION_BRANCH="ftk-v1-gpu-acceleration"
SOURCE_BRANCH="no_ndarray"
BASE_BRANCH="master"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Usage
usage() {
    echo "Usage: $0 [step]"
    echo ""
    echo "Steps:"
    echo "  1 - Create integration branch and backup"
    echo "  2 - Cherry-pick infrastructure improvements"
    echo "  3 - Integrate GPU-specific code"
    echo "  4 - Update build system"
    echo "  5 - Testing and validation"
    echo "  all - Run all steps (use with caution!)"
    exit 1
}

# Step 1: Create integration branch
step1_create_branch() {
    echo -e "${BLUE}=== Step 1: Creating Integration Branch ===${NC}"
    echo ""

    # Check if on master
    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "$BASE_BRANCH" ]; then
        echo -e "${YELLOW}Warning: Not on $BASE_BRANCH branch (currently on $current_branch)${NC}"
        read -p "Switch to $BASE_BRANCH? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            git checkout "$BASE_BRANCH"
        else
            echo "Aborting."
            exit 1
        fi
    fi

    # Ensure working directory is clean
    if ! git diff-index --quiet HEAD --; then
        echo -e "${RED}Error: Working directory is not clean${NC}"
        echo "Please commit or stash your changes first."
        exit 1
    fi

    # Create backup tag
    backup_tag="backup-master-$(date +%Y%m%d-%H%M%S)"
    git tag "$backup_tag"
    echo -e "${GREEN}Created backup tag: $backup_tag${NC}"

    # Create integration branch
    if git rev-parse --verify "$INTEGRATION_BRANCH" >/dev/null 2>&1; then
        echo -e "${YELLOW}Integration branch $INTEGRATION_BRANCH already exists${NC}"
        read -p "Delete and recreate? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            git branch -D "$INTEGRATION_BRANCH"
        else
            echo "Using existing branch"
            git checkout "$INTEGRATION_BRANCH"
            return
        fi
    fi

    git checkout -b "$INTEGRATION_BRANCH"
    echo -e "${GREEN}Created integration branch: $INTEGRATION_BRANCH${NC}"
    echo ""
    echo "Next steps:"
    echo "  ./scripts/integrate_no_ndarray.sh 2"
}

# Step 2: Cherry-pick infrastructure improvements
step2_infrastructure() {
    echo -e "${BLUE}=== Step 2: Cherry-picking Infrastructure Improvements ===${NC}"
    echo ""

    # Ensure on integration branch
    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "$INTEGRATION_BRANCH" ]; then
        echo -e "${RED}Error: Not on $INTEGRATION_BRANCH branch${NC}"
        exit 1
    fi

    # Key infrastructure files to integrate
    echo "Integrating key infrastructure files from $SOURCE_BRANCH:"
    echo ""

    # 1. kd_lite - GPU-compatible KD-tree
    echo "1. Integrating kd_lite (GPU-compatible spatial indexing)..."
    if git show "$SOURCE_BRANCH:include/ftk/basic/kd_lite.hh" > /dev/null 2>&1; then
        git checkout "$SOURCE_BRANCH" -- include/ftk/basic/kd_lite.hh
        git add include/ftk/basic/kd_lite.hh
        echo -e "   ${GREEN}✓ kd_lite.hh added${NC}"
    else
        echo -e "   ${YELLOW}⚠ kd_lite.hh not found in $SOURCE_BRANCH${NC}"
    fi

    # 2. MPAS particle structure
    echo "2. Integrating MPAS particle structure..."
    if git show "$SOURCE_BRANCH:include/ftk/features/mpas_particle.hh" > /dev/null 2>&1; then
        git checkout "$SOURCE_BRANCH" -- include/ftk/features/mpas_particle.hh
        git add include/ftk/features/mpas_particle.hh
        echo -e "   ${GREEN}✓ mpas_particle.hh added${NC}"
    else
        echo -e "   ${YELLOW}⚠ mpas_particle.hh not found${NC}"
    fi

    # 3. Error handling improvements
    echo "3. Integrating error handling improvements..."
    if git show "$SOURCE_BRANCH:include/ftk/error.hh" > /dev/null 2>&1; then
        # Check for conflicts first
        git show "$BASE_BRANCH:include/ftk/error.hh" > /tmp/base_error.hh 2>/dev/null || touch /tmp/base_error.hh
        git show "$SOURCE_BRANCH:include/ftk/error.hh" > /tmp/source_error.hh

        if diff -q /tmp/base_error.hh /tmp/source_error.hh > /dev/null 2>&1; then
            echo -e "   ${YELLOW}⚠ No changes in error.hh${NC}"
        else
            git checkout "$SOURCE_BRANCH" -- include/ftk/error.hh
            git add include/ftk/error.hh
            echo -e "   ${GREEN}✓ error.hh updated${NC}"
        fi

        rm -f /tmp/base_error.hh /tmp/source_error.hh
    fi

    # 4. Updated JSON library
    echo "4. Checking JSON library updates..."
    # Usually external libraries should be handled carefully
    echo -e "   ${YELLOW}⚠ Skipping automatic merge of external/json.hh (needs manual review)${NC}"

    # Commit infrastructure changes
    if git diff --cached --quiet; then
        echo ""
        echo -e "${YELLOW}No files staged for commit${NC}"
    else
        git commit -m "Integrate GPU infrastructure from no_ndarray branch

- Add kd_lite: GPU-compatible KD-tree implementation
- Add mpas_particle_t: Spherical coordinate particle structure
- Update error handling

Source: $SOURCE_BRANCH branch
Part of FTK v1.0 GPU acceleration initiative"

        echo ""
        echo -e "${GREEN}Infrastructure changes committed${NC}"
    fi

    echo ""
    echo "Next steps:"
    echo "  ./scripts/integrate_no_ndarray.sh 3"
}

# Step 3: Integrate GPU-specific code
step3_gpu_code() {
    echo -e "${BLUE}=== Step 3: Integrating GPU-Specific Code ===${NC}"
    echo ""

    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "$INTEGRATION_BRANCH" ]; then
        echo -e "${RED}Error: Not on $INTEGRATION_BRANCH branch${NC}"
        exit 1
    fi

    # List GPU files from source branch
    echo "GPU files in $SOURCE_BRANCH:"
    git ls-tree -r --name-only "$SOURCE_BRANCH" | grep -E '\.(cu|cuh)$' || echo "  None found"
    echo ""

    # MPAS ocean particle tracker CUDA
    echo "Integrating MPAS ocean GPU tracker..."
    if git show "$SOURCE_BRANCH:include/ftk/filters/mpas_ocean_particle_tracker.cuh" > /dev/null 2>&1; then
        # Create directory if needed
        mkdir -p include/ftk/filters

        git checkout "$SOURCE_BRANCH" -- include/ftk/filters/mpas_ocean_particle_tracker.cuh
        git add include/ftk/filters/mpas_ocean_particle_tracker.cuh
        echo -e "   ${GREEN}✓ mpas_ocean_particle_tracker.cuh added${NC}"
    else
        echo -e "   ${YELLOW}⚠ mpas_ocean_particle_tracker.cuh not found${NC}"
    fi

    # MPAS particle tracer header
    if git show "$SOURCE_BRANCH:include/ftk/filters/particle_tracer_mpas_ocean.hh" > /dev/null 2>&1; then
        git checkout "$SOURCE_BRANCH" -- include/ftk/filters/particle_tracer_mpas_ocean.hh
        git add include/ftk/filters/particle_tracer_mpas_ocean.hh
        echo -e "   ${GREEN}✓ particle_tracer_mpas_ocean.hh added${NC}"
    fi

    # Commit GPU code
    if git diff --cached --quiet; then
        echo ""
        echo -e "${YELLOW}No GPU files staged${NC}"
    else
        git commit -m "Add GPU-accelerated MPAS ocean particle tracker

- CUDA implementation of particle advection
- Spherical coordinate RK4 integration on GPU
- Spatial indexing with kd_lite

Source: $SOURCE_BRANCH branch
Performance: Expected 10-50x speedup for large particle counts"

        echo ""
        echo -e "${GREEN}GPU code committed${NC}"
    fi

    echo ""
    echo "Next steps:"
    echo "  ./scripts/integrate_no_ndarray.sh 4"
}

# Step 4: Update build system
step4_build_system() {
    echo -e "${BLUE}=== Step 4: Updating Build System ===${NC}"
    echo ""

    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "$INTEGRATION_BRANCH" ]; then
        echo -e "${RED}Error: Not on $INTEGRATION_BRANCH branch${NC}"
        exit 1
    fi

    echo "Creating CMake option for GPU-native mode..."

    # Backup CMakeLists.txt
    cp CMakeLists.txt CMakeLists.txt.backup

    # Add GPU-native option (manual edit recommended)
    cat << 'EOF' > /tmp/ftk_cmake_addition.txt

# GPU-native mode (lightweight data structures, no ndarray overhead)
option(FTK_GPU_NATIVE "Use GPU-native data structures (no ndarray)" OFF)

if(FTK_GPU_NATIVE)
  message(STATUS "FTK: GPU-native mode enabled")
  add_definitions(-DFTK_GPU_NATIVE)
endif()
EOF

    echo ""
    echo "CMake additions prepared in /tmp/ftk_cmake_addition.txt"
    echo ""
    echo -e "${YELLOW}⚠ Manual CMakeLists.txt edit required:${NC}"
    echo "  1. Open CMakeLists.txt in your editor"
    echo "  2. Add the content from /tmp/ftk_cmake_addition.txt after line ~13 (after FTK_BUILD_XGC_UTILS option)"
    echo "  3. Update config.hh.in to add FTK_GPU_NATIVE define"
    echo ""
    read -p "Press Enter when done editing..."

    # Update config.hh.in
    if [ -f include/ftk/config.hh.in ]; then
        echo "Updating config.hh.in..."
        if ! grep -q "FTK_GPU_NATIVE" include/ftk/config.hh.in; then
            echo "#cmakedefine01 FTK_GPU_NATIVE" >> include/ftk/config.hh.in
            echo -e "   ${GREEN}✓ config.hh.in updated${NC}"
        else
            echo -e "   ${YELLOW}⚠ FTK_GPU_NATIVE already in config.hh.in${NC}"
        fi
    fi

    # Commit build system changes
    git add CMakeLists.txt include/ftk/config.hh.in
    git commit -m "Add GPU-native build mode option

- Add FTK_GPU_NATIVE CMake option
- Configure GPU-native data structures
- Prepare for no_ndarray architecture

This allows toggling between:
  - Legacy: ftk::ndarray-based (better for CPU, Python, VTK)
  - GPU-native: Lightweight structures (better for GPU performance)" || echo "No changes to commit"

    echo ""
    echo -e "${GREEN}Build system updated${NC}"
    echo ""
    echo "Next steps:"
    echo "  ./scripts/integrate_no_ndarray.sh 5"
}

# Step 5: Testing
step5_testing() {
    echo -e "${BLUE}=== Step 5: Testing and Validation ===${NC}"
    echo ""

    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "$INTEGRATION_BRANCH" ]; then
        echo -e "${RED}Error: Not on $INTEGRATION_BRANCH branch${NC}"
        exit 1
    fi

    # Build test
    echo "Building FTK with GPU-native mode..."
    if [ -d build ]; then
        echo -e "${YELLOW}Removing old build directory...${NC}"
        rm -rf build
    fi

    mkdir -p build && cd build

    echo ""
    echo "Configuring with GPU-native mode enabled..."
    cmake .. \
        -DFTK_GPU_NATIVE=ON \
        -DFTK_USE_CUDA=ON \
        -DFTK_BUILD_TESTS=ON \
        -DCMAKE_BUILD_TYPE=Release

    echo ""
    echo "Building..."
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ Build successful!${NC}"
    else
        echo ""
        echo -e "${RED}✗ Build failed${NC}"
        cd ..
        exit 1
    fi

    # Run tests
    if [ -f bin/ftk ]; then
        echo ""
        echo "Testing FTK binary..."
        ./bin/ftk --help
    fi

    cd ..

    echo ""
    echo -e "${GREEN}=== Integration Complete! ===${NC}"
    echo ""
    echo "Summary:"
    echo "  Branch: $INTEGRATION_BRANCH"
    git log --oneline master..$INTEGRATION_BRANCH
    echo ""
    echo "Next steps:"
    echo "  1. Run comprehensive tests: cd build && ctest"
    echo "  2. Test with real datasets (MPAS, XGC)"
    echo "  3. Benchmark GPU vs CPU performance"
    echo "  4. If satisfied: git checkout master && git merge $INTEGRATION_BRANCH"
}

# Main
STEP=${1:-}

case "$STEP" in
    1)
        step1_create_branch
        ;;
    2)
        step2_infrastructure
        ;;
    3)
        step3_gpu_code
        ;;
    4)
        step4_build_system
        ;;
    5)
        step5_testing
        ;;
    all)
        echo -e "${YELLOW}Warning: Running all steps automatically${NC}"
        read -p "Are you sure? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            step1_create_branch
            step2_infrastructure
            step3_gpu_code
            step4_build_system
            step5_testing
        fi
        ;;
    *)
        usage
        ;;
esac
