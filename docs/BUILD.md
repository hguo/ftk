# FTK Build System Documentation

This document provides comprehensive guidance for building FTK (Feature Tracking Kit) with various configurations and dependencies.

## Table of Contents
- [Quick Start](#quick-start)
- [Build Requirements](#build-requirements)
- [CMake Options](#cmake-options)
- [Minimal Build](#minimal-build)
- [Feature Builds](#feature-builds)
- [Common Build Configurations](#common-build-configurations)
- [Troubleshooting](#troubleshooting)

## Quick Start

### Minimal Build (Core Features Only)

```bash
git clone https://github.com/hguo/ftk.git
cd ftk
mkdir build && cd build
cmake ..
make -j8
make install
```

### Recommended Build (Most Features)

```bash
cmake .. \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_BUILD_TESTS=ON
make -j8
ctest
make install
```

## Build Requirements

### Required Dependencies
- **CMake** >= 3.10
- **C++17 compliant compiler** (GCC 7+, Clang 5+, MSVC 2017+)
- **ndarray library** (required, external dependency)
- **yaml-cpp** (for configuration file parsing)

### Optional Dependencies
See [CMake Options](#cmake-options) below for the full list of optional dependencies.

## CMake Options

FTK uses a flexible option system where each dependency can be set to:
- **TRUE**: Require the dependency (fail if not found)
- **AUTO**: Use if available, skip if not found
- **FALSE**: Disable the dependency (default for most)

### Core Build Options

| Option | Description | Default |
|--------|-------------|---------|
| `FTK_BUILD_TESTS` | Build test suite | OFF |
| `FTK_BUILD_PARAVIEW` | Build ParaView plugins | OFF |
| `FTK_BUILD_PYFTK` | Build Python bindings | OFF |
| `FTK_BUILD_XGC_UTILS` | Build XGC utilities | OFF |

### Optional Feature Dependencies

| Option | Description | Status | Use Case |
|--------|-------------|--------|----------|
| `FTK_USE_Boost` | Boost C++ Libraries | Not actively used | Legacy support |
| `FTK_USE_CGAL` | Computational Geometry Algorithms Library | Experimental | Advanced mesh operations |
| `FTK_USE_CUDA` | NVIDIA CUDA support | Experimental | GPU acceleration |
| `FTK_USE_DECAF` | Decaf data flow framework | Experimental | In-situ workflows |
| `FTK_USE_GSL` | GNU Scientific Library | Optional | Numerical algorithms |
| `FTK_USE_GMP` | GNU Multiple Precision Arithmetic | Optional | Exact arithmetic (required by CGAL) |
| `FTK_USE_Kokkos` | Kokkos performance portability | Experimental | Multi-architecture support |
| `FTK_USE_LevelDB` | LevelDB key-value store | Optional | Feature persistence |
| `FTK_USE_METIS` | METIS graph partitioning | Optional | Domain decomposition |
| `FTK_USE_MPFR` | Multiple Precision Floating-Point | Optional | High-precision arithmetic (required by CGAL) |
| `FTK_USE_MPI` | Message Passing Interface | **Recommended** | Distributed computing |
| `FTK_USE_MPSolve` | Multiprecision Polynomial Solver | Optional | Polynomial root finding |
| `FTK_USE_OpenMP` | OpenMP threading | **Recommended** | Shared-memory parallelism |
| `FTK_USE_Qt5` | Qt5 GUI framework | Optional | Visualization tools |
| `FTK_USE_RocksDB` | RocksDB key-value store | Optional | Feature persistence |
| `FTK_USE_HIPSYCL` | hipSYCL (SYCL for AMD/NVIDIA) | Experimental | GPU acceleration |
| `FTK_USE_SYCL` | SYCL standard | Experimental | Cross-platform GPU support |
| `FTK_USE_TBB` | Intel Thread Building Blocks | **Recommended** | Advanced threading |
| `FTK_USE_VTK` | Visualization Toolkit | **Recommended** | Mesh I/O and visualization |

## Minimal Build

The absolute minimum build requires only core dependencies:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local
make -j$(nproc)
make install
```

**Capabilities:**
- ✅ Core tracker algorithms
- ✅ Regular grid tracking
- ✅ Basic I/O (text, binary)
- ✅ Single-threaded execution
- ❌ No MPI support
- ❌ No VTK I/O
- ❌ No GPU acceleration

## Feature Builds

### Standard Scientific Computing Build

For most scientific computing workflows:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_TBB=AUTO \
  -DFTK_BUILD_TESTS=ON
```

**Provides:**
- Distributed computing with MPI
- Shared-memory parallelism with OpenMP
- VTK mesh I/O and visualization
- Optional TBB for advanced threading
- Test suite for verification

### High-Performance Computing (HPC) Build

For large-scale HPC systems:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_TBB=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_METIS=TRUE \
  -DCMAKE_INSTALL_PREFIX=$HOME/ftk
```

**Provides:**
- Full MPI+OpenMP hybrid parallelism
- TBB for fine-grained parallelism
- METIS for load balancing
- VTK for data I/O

### GPU-Accelerated Build

For systems with NVIDIA GPUs:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DFTK_USE_CUDA=TRUE \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DCMAKE_CUDA_ARCHITECTURES="70;80" \
  -DCMAKE_CUDA_COMPILER=nvcc
```

**Provides:**
- CUDA acceleration for critical point detection
- GPU-based feature tracking
- MPI for multi-GPU scaling

### XGC Fusion Simulation Build

For XGC plasma simulation data:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_TBB=TRUE \
  -DFTK_BUILD_XGC_UTILS=ON \
  -DFTK_USE_ADIOS2=TRUE \
  -DFTK_USE_HDF5=TRUE
```

**Provides:**
- XGC mesh support
- ADIOS2 and HDF5 I/O
- Specialized XGC trackers

### ParaView Plugin Build

For interactive visualization:

```bash
cmake .. \
  -DFTK_BUILD_PARAVIEW=ON \
  -DParaView_DIR=/path/to/paraview/lib/cmake/paraview-5.10 \
  -DCMAKE_BUILD_TYPE=Release
```

**Provides:**
- ParaView plugins for FTK algorithms
- Interactive feature tracking
- Real-time visualization

### Python Bindings Build

For Python users:

```bash
cmake .. \
  -DFTK_BUILD_PYFTK=ON \
  -DPYTHON_EXECUTABLE=$(which python3) \
  -DCMAKE_BUILD_TYPE=Release
```

**Provides:**
- Python module `pyftk`
- NumPy integration
- Jupyter notebook support

## Common Build Configurations

### Development Build

For FTK developers:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-g -O0 -Wall -Wextra" \
  -DFTK_BUILD_TESTS=ON \
  -DFTK_USE_VTK=TRUE \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

### AddressSanitizer Build

For memory leak detection:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address" \
  -DFTK_BUILD_TESTS=ON
```

### ThreadSanitizer Build

For thread-safety testing:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_BUILD_TESTS=ON
```

## Platform-Specific Notes

### macOS

```bash
# Install dependencies with Homebrew
brew install cmake vtk tbb open-mpi libomp

cmake .. \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_TBB=TRUE \
  -DOpenMP_ROOT=$(brew --prefix libomp)
```

### Linux (Ubuntu/Debian)

```bash
# Install dependencies
sudo apt-get install cmake build-essential \
  libvtk9-dev libtbb-dev libopenmpi-dev \
  libyaml-cpp-dev

cmake .. \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_TBB=TRUE
```

### Linux (RHEL/CentOS)

```bash
# Install dependencies
sudo yum install cmake3 gcc-c++ \
  vtk-devel tbb-devel openmpi-devel \
  yaml-cpp-devel

module load mpi
cmake3 .. \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE
```

### HPC Systems (Example: NERSC Cori)

```bash
module load cmake vtk tbb

cmake .. \
  -DCMAKE_CXX_COMPILER=CC \
  -DFTK_USE_MPI=TRUE \
  -DFTK_USE_OpenMP=TRUE \
  -DFTK_USE_VTK=TRUE \
  -DFTK_USE_TBB=TRUE \
  -DCMAKE_INSTALL_PREFIX=$SCRATCH/ftk
```

## Troubleshooting

### CMake Cannot Find Dependencies

**Problem**: `Could not find package VTK`

**Solutions**:
1. Set package hints:
   ```bash
   cmake .. -DVTK_DIR=/path/to/vtk/lib/cmake/vtk-9.0
   ```

2. Use AUTO mode to skip if not found:
   ```bash
   cmake .. -DFTK_USE_VTK=AUTO
   ```

3. Install the dependency:
   ```bash
   # Ubuntu
   sudo apt-get install libvtk9-dev

   # macOS
   brew install vtk
   ```

### Compiler Version Too Old

**Problem**: `error: #error This file requires compiler and library support for the ISO C++ 2017 standard`

**Solution**: Use a newer compiler:
```bash
cmake .. -DCMAKE_CXX_COMPILER=g++-9
```

Or on HPC systems:
```bash
module load gcc/9.3.0
```

### MPI Compilation Errors

**Problem**: Cannot find MPI headers

**Solution**: Use MPI compiler wrapper:
```bash
cmake .. -DCMAKE_CXX_COMPILER=mpicxx
```

Or specify MPI paths:
```bash
cmake .. \
  -DMPI_CXX_COMPILER=/usr/bin/mpicxx \
  -DFTK_USE_MPI=TRUE
```

### CUDA Build Failures

**Problem**: `nvcc fatal : Unsupported gpu architecture 'compute_80'`

**Solution**: Specify compatible architectures:
```bash
cmake .. \
  -DFTK_USE_CUDA=TRUE \
  -DCMAKE_CUDA_ARCHITECTURES="60;70;75"  # Adjust for your GPU
```

### Out of Memory During Compilation

**Problem**: Compiler crashes or system freezes

**Solution**: Reduce parallel jobs:
```bash
make -j4  # Instead of -j$(nproc)
```

Or build specific targets:
```bash
make libftk  # Build library only
make ftk     # Build executable only
```

### Test Failures

**Problem**: `ctest` reports failures

**Solution**: Run specific test with verbose output:
```bash
ctest -R test_critical_point_tracking_double_gyre -VV
```

Or run test binary directly:
```bash
./bin/test_critical_point_tracking_double_gyre
```

## Build Performance Tips

### Faster Builds

1. **Use Ninja instead of Make**:
   ```bash
   cmake .. -G Ninja
   ninja -j$(nproc)
   ```

2. **Use ccache** for incremental builds:
   ```bash
   cmake .. -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
   ```

3. **Build only what you need**:
   ```bash
   make libftk ftk  # Skip tests and examples
   ```

### Smaller Binaries

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=MinSizeRel \
  -DCMAKE_CXX_FLAGS="-Os -flto"
```

## Environment Variables

Useful environment variables for controlling the build:

| Variable | Description |
|----------|-------------|
| `CMAKE_PREFIX_PATH` | Search paths for dependencies |
| `CMAKE_INSTALL_PREFIX` | Installation directory |
| `CMAKE_BUILD_TYPE` | Debug, Release, RelWithDebInfo, MinSizeRel |
| `CMAKE_CXX_COMPILER` | C++ compiler to use |
| `CMAKE_CXX_FLAGS` | Additional compiler flags |
| `OMP_NUM_THREADS` | OpenMP thread count |
| `CUDA_VISIBLE_DEVICES` | GPU selection for CUDA builds |

## Verifying the Build

After building, verify functionality:

```bash
# Run tests
ctest --output-on-failure

# Check binary
./bin/ftk --help

# Verify library
nm lib/libftk.a | grep critical_point_tracker

# Python bindings (if built)
python3 -c "import pyftk; print(pyftk.__version__)"
```

## Next Steps

After a successful build:
1. Run `make install` to install to `CMAKE_INSTALL_PREFIX`
2. Read [THREAD_SAFETY.md](THREAD_SAFETY.md) for parallel execution guidance
3. See examples in `tests/` directory
4. Generate API documentation: `doxygen Doxyfile`

## Getting Help

- **Build issues**: https://github.com/hguo/ftk/issues
- **Documentation**: https://github.com/hguo/ftk/wiki
- **Mailing list**: (if available)

---

**Last Updated**: February 2026
**FTK Version**: Post ndarray migration (no_ndarray branch)
