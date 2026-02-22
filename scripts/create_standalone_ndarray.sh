#!/bin/bash
# create_standalone_ndarray.sh - Extract FTK's old ndarray into standalone library
#
# This script creates a temporary standalone ndarray library from FTK's history
# until an official hguo/ndarray repository is created.
#
# Usage: ./scripts/create_standalone_ndarray.sh [output_dir]

set -e

OUTPUT_DIR=${1:-"$HOME/ndarray-standalone"}
BACKUP_TAG="backup-master-before-no_ndarray-merge-20260115-154320"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}=== Creating Standalone ndarray Library ===${NC}"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if we're in FTK repo
if [ ! -f "CMakeLists.txt" ] || ! grep -q "project (ftk)" CMakeLists.txt; then
    echo -e "${RED}Error: Must run from FTK root directory${NC}"
    exit 1
fi

# Save current state
CURRENT_BRANCH=$(git branch --show-current)
echo "Current branch: $CURRENT_BRANCH"

# Check if backup tag exists
if ! git rev-parse --verify "$BACKUP_TAG" >/dev/null 2>&1; then
    echo -e "${RED}Error: Backup tag $BACKUP_TAG not found${NC}"
    echo "This script needs the pre-merge FTK version with internal ndarray"
    exit 1
fi

# Check if output directory exists
if [ -d "$OUTPUT_DIR" ]; then
    echo -e "${YELLOW}Warning: $OUTPUT_DIR already exists${NC}"
    read -p "Delete and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$OUTPUT_DIR"
    else
        echo "Aborting."
        exit 1
    fi
fi

mkdir -p "$OUTPUT_DIR"
echo -e "${GREEN}Created $OUTPUT_DIR${NC}"

# Extract ndarray from old FTK
echo ""
echo "Extracting ndarray from FTK history..."
git checkout "$BACKUP_TAG" -- include/ftk/ndarray 2>/dev/null || {
    echo -e "${RED}Error: Could not extract ndarray directory${NC}"
    exit 1
}

# Copy to output directory
cp -r include/ftk/ndarray "$OUTPUT_DIR/"
echo -e "${GREEN}✓ Copied ndarray files${NC}"

# Also copy the main ndarray.hh if it exists
if git show "$BACKUP_TAG:include/ftk/ndarray.hh" > "$OUTPUT_DIR/ndarray/ndarray.hh" 2>/dev/null; then
    echo -e "${GREEN}✓ Copied ndarray.hh${NC}"
fi

# Clean up git changes
git reset HEAD include/ftk/ndarray 2>/dev/null || true
git checkout -- include/ftk/ndarray 2>/dev/null || true

# Create CMakeLists.txt for standalone library
cat > "$OUTPUT_DIR/CMakeLists.txt" << 'EOF'
cmake_minimum_required(VERSION 3.10)
project(ndarray VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options for I/O support (from original FTK)
option(NDARRAY_USE_VTK "Enable VTK support" OFF)
option(NDARRAY_USE_NETCDF "Enable NetCDF support" OFF)
option(NDARRAY_USE_HDF5 "Enable HDF5 support" OFF)
option(NDARRAY_USE_ADIOS2 "Enable ADIOS2 support" OFF)

# Header-only interface library
add_library(ndarray INTERFACE)

target_include_directories(ndarray INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Dependencies
if(NDARRAY_USE_VTK)
    find_package(VTK REQUIRED)
    target_link_libraries(ndarray INTERFACE ${VTK_LIBRARIES})
    target_compile_definitions(ndarray INTERFACE NDARRAY_HAVE_VTK)
endif()

if(NDARRAY_USE_NETCDF)
    find_package(NetCDF REQUIRED)
    target_link_libraries(ndarray INTERFACE NetCDF::NetCDF)
    target_compile_definitions(ndarray INTERFACE NDARRAY_HAVE_NETCDF)
endif()

if(NDARRAY_USE_HDF5)
    find_package(HDF5 REQUIRED)
    target_link_libraries(ndarray INTERFACE ${HDF5_LIBRARIES})
    target_include_directories(ndarray INTERFACE ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(ndarray INTERFACE NDARRAY_HAVE_HDF5)
endif()

if(NDARRAY_USE_ADIOS2)
    find_package(ADIOS2 REQUIRED)
    target_link_libraries(ndarray INTERFACE adios2::adios2)
    target_compile_definitions(ndarray INTERFACE NDARRAY_HAVE_ADIOS2)
endif()

# Installation
install(TARGETS ndarray EXPORT ndarrayTargets)
install(DIRECTORY ndarray/ DESTINATION include/ndarray
        FILES_MATCHING PATTERN "*.hh")

# CMake package config
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ndarrayConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(EXPORT ndarrayTargets
    FILE ndarrayTargets.cmake
    NAMESPACE ndarray::
    DESTINATION lib/cmake/ndarray
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ndarrayConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ndarrayConfig.cmake"
    INSTALL_DESTINATION lib/cmake/ndarray
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ndarrayConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ndarrayConfigVersion.cmake"
    DESTINATION lib/cmake/ndarray
)

# Print configuration summary
message(STATUS "")
message(STATUS "ndarray configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  VTK support: ${NDARRAY_USE_VTK}")
message(STATUS "  NetCDF support: ${NDARRAY_USE_NETCDF}")
message(STATUS "  HDF5 support: ${NDARRAY_USE_HDF5}")
message(STATUS "  ADIOS2 support: ${NDARRAY_USE_ADIOS2}")
message(STATUS "")
EOF

echo -e "${GREEN}✓ Created CMakeLists.txt${NC}"

# Create cmake directory and config
mkdir -p "$OUTPUT_DIR/cmake"
cat > "$OUTPUT_DIR/cmake/ndarrayConfig.cmake.in" << 'EOF'
@PACKAGE_INIT@

include("${CMAKE_CURRENT_LIST_DIR}/ndarrayTargets.cmake")

check_required_components(ndarray)
EOF

echo -e "${GREEN}✓ Created CMake package config${NC}"

# Create README
cat > "$OUTPUT_DIR/README.md" << 'EOF'
# ndarray - N-Dimensional Array Library

**Extracted from FTK (Feature Tracking Kit)**
**Version**: 0.1.0 (temporary standalone)
**Purpose**: Provide ndarray functionality for FTK v1.0+

## Overview

This is a temporary standalone version of FTK's ndarray implementation,
extracted to satisfy FTK's external ndarray dependency.

**Status**: This is a temporary solution until an official hguo/ndarray
repository is created and published.

## Features

- N-dimensional array container (row-major order)
- I/O support: NetCDF, HDF5, ADIOS2, VTK
- Gradient computation
- Synthetic data generation
- Stream processing

## Installation

### Basic Installation (no I/O dependencies)

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
```

### With I/O Support

```bash
cmake .. \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DNDARRAY_USE_NETCDF=ON \
  -DNDARRAY_USE_HDF5=ON \
  -DNDARRAY_USE_ADIOS2=ON \
  -DNDARRAY_USE_VTK=ON
make
sudo make install
```

## Usage with FTK

After installing ndarray, build FTK:

```bash
cd /path/to/ftk
mkdir build && cd build
cmake .. -Dndarray_DIR=/usr/local/lib/cmake/ndarray
make
```

## Documentation

See original FTK documentation:
- [FTK ndarray docs](https://github.com/hguo/ftk/blob/master/docs/ndarray.md)

## License

MIT License (inherited from FTK)

## TODO

- [ ] Create official hguo/ndarray repository
- [ ] Improve documentation
- [ ] Add examples
- [ ] Add tests
- [ ] Optimize for modern C++
- [ ] Add GPU support

## Credits

Extracted from FTK: https://github.com/hguo/ftk
Original author: Hanqi Guo

EOF

echo -e "${GREEN}✓ Created README.md${NC}"

# Create .gitignore
cat > "$OUTPUT_DIR/.gitignore" << 'EOF'
build/
*.swp
*.swo
*~
.DS_Store
EOF

# Print summary
echo ""
echo -e "${GREEN}=== Standalone ndarray Created! ===${NC}"
echo ""
echo "Location: $OUTPUT_DIR"
echo ""
echo "Files created:"
find "$OUTPUT_DIR" -type f -name "*.hh" | wc -l | xargs echo "  - Header files:"
echo "  - CMakeLists.txt"
echo "  - cmake/ndarrayConfig.cmake.in"
echo "  - README.md"
echo ""
echo "Next steps:"
echo ""
echo "1. Build and install ndarray:"
echo "   cd $OUTPUT_DIR"
echo "   mkdir build && cd build"
echo "   cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local"
echo "   make"
echo "   sudo make install"
echo ""
echo "2. Build FTK with ndarray:"
echo "   cd $(pwd)"
echo "   mkdir build && cd build"
echo "   cmake .. -Dndarray_DIR=/usr/local/lib/cmake/ndarray"
echo "   make"
echo ""
echo "3. (Optional) Create git repository:"
echo "   cd $OUTPUT_DIR"
echo "   git init"
echo "   git add ."
echo "   git commit -m 'Initial commit: Standalone ndarray from FTK'"
echo "   # Then push to GitHub as hguo/ndarray"
echo ""
