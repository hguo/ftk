# FTK Installation

## Spack install

One can install FTK through [spack](https://spack.io/) (we recommend to use the latest version of spack on the develop branch):  

```bash
$ spack install ftk
$ spack load ftk
$ ftk
Usage:
  ftk [OPTION...]
  -f, --feature arg             Feature type: critical_point, isosurface,
  .....
```

To build FTK's advanced features, one can build and install FTK with various optional dependencies including ADIOS2, CUDA, MPI, NetCDF, and VTK; use the `master` branch for the latest features, for example:

```bash
$ spack install ftk+vtk+netcdf+mpi+cuda+adios2@master
```

See `spack info ftk` for more details.

## Pip install (Python bindings only)

One can install PyFTK through PyPI: 

```bash
$ pip3 install pyftk
```

## Build and install FTK from source

**Important note:** If you change CMake build options, e.g. `FTK_USE_HDF5`, we recommend to reconfigure with CMake in an empty build directory.  If you have previously installed FTK to a system directory e.g. `/usr/local` , it is important to completely delete all FTK headers and binaries; otherwise your compiler might find FTK headers in default search directories such as `/usr/local/include/ftk/config.hh`, causing unexpected errors. 

### Dependencies and CMake options

Mandatory dependencies

* CMake (3.10 minimum)

Optional dependencies

* ADIOS2 (>=2.7.0): Use `-DFTK_USE_ADIOS2=ON` with cmake; we suggest to use MPI when building with ADIOS2
* CUDA (>=10.1): Use `-DFTK_USE_CUDA=ON` and specify `-DCMAKE_CUDA_COMPILER`
* HDF5: Use `-DFTK_USE_HDF5=ON`
* GMP (strongly recommended): will automatically find GMP unless `-DFTK_USE_GMP=OFF`
* MPI (mpich 3.2 and higher versions recommended): use `CC=mpicc CXX=mpicxx` and `-DFTK_USE_MPI=ON`, see more details below
* MPSolve: Use `-DFTK_USE_MPSOLVE=ON`
* NetCDF-C: Use `-DFTK_USE_NETCDF=ON` and specify `-DNETCDF_DIR`
* ParaView (5.8.0 and higher versions recommended): use `-DFTK_BUILD_PARAVIEW` and specify `-DParaView_DIR` in cmake options; see more details below
* Python3: use `setup.py` or `-DFTK_BUILD_PYFTK=ON` in make options; see more details below
* VTK (9.0.1 and higher versions recommended): use `-DFTK_BUILD_VTK` and specify `-DVTK_DIR`/`-DCMAKE_PREFIX_PATH` in cmake options; see more details below

###  Building without external dependencies

```bash
$ cd $FTK_SOURCE_DIR/build
$ cmake .. && make && make install
```

### Building with MPI

You may use MPI to accelerate feature tracking with both distributed-parallelism.  To build FTK with MPI, you need to use MPI C/C++ compilers: 

```bash
$ CC=mpicc CXX=mpicxx cmake -DFTK_USE_MPI=ON
```

### Building ParaView plugins

FTK provides ParaView plugins to allow users track critical points (maxima, minima, and saddles) in scalar field data.  In order to build the plugins, one needs to build ParaView from source, for example: 

```bash
$ tar zxf ParaView-v5.8.0.tar.gz
$ cd ParaView-v5.8.0
$ mkdir build && cd build
$ cmake ..   # You may need to specify paths to Qt and other dependencies needed by ParaView
$ make
```

Assuming `$YOUR_ParaView_Build` is the above build directory, one can further build ParaView plugins as follows

```bash
$ git clone https://github.com/hguo/ftk $FTK_SOURCE_DIR
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_PARAVIEW=ON -DParaView_DIR=$YOUR_ParaView_Build
$ make
```

If built successfully, you will see the plugins binary as `lib/paraview-5.8/plugins/FTK/FTK.so`.  Open the "Plugin Manager" in ParaView GUI, and load this binary with "Load New..." button, and then select and load FTK in the list.  To check if ParaView plugins are correctly built by reproducing the results in the above figure, use "Sources-->FTK-->SpiralWoven2DSource", "Filters-->FTK-->CriticalPointTracker2D",followed by the "Tube" filter in ParaView.

### Building PyFTK

FTK Python bindings is based on [pybind11](https://github.com/pybind/pybind11), which is included as a submodule under the `python` subdirectory; make sure the submodule is checked out:

```bash
$ git submodule update --init --recursive
```

If you downloaded a release from tarball archives on GitHub, be sure to manually download `pybind11` to the designated directory before building.  

You may build PyFTK with `setuptools` or CMake.  Notice that CMake is required to build PyFTK.  We do not recommend to use Python 2.  Advanced build options is currently not possible to configure with `setuptools`.  

Build PyFTK with `setuptools`:

```bash
$ cd $FTK_SOURCE_DIR
$ python setup.py install
```

Build PyFTK with CMake:

```bash
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_PYFTK=ON
$ make
```

The output PyFTK binary will be in the `lib` directory.
