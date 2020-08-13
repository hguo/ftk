# FTK: The Feature Tracking Kit

[![Build Status](https://travis-ci.org/hguo/ftk.svg?branch=master)](https://travis-ci.org/hguo/ftk)

FTK is a library that provides building blocks for feature tracking algorithms in scientific datasets.  

![](docs/images/critical_point_tracking_2d_paraview.png)

## Major components in FTK

* Hypermesh: data structures for high-dimensional meshes and mesh elements including *n*-simplices, *n*-cubes, and *n*-prisms; utilities to generalize given 2D/3D structured/unstructured meshes into 3D/4D spacetime meshes

![](./docs/images/regular_simplex_subdivision.svg)

* Numeric: root-find algorithms for inverse interpolations and parallel vector operators in *n*-simplices, *n*-cubes, and simplex-prisms; lightweight linear algebra utilities to support root-finding
* CCL: connected component labeling algorithm for building feature tracking algorithms
* Geometry: utilities to transform connect components to geometry for visualization and analysis
* Tracking graph: data structures to record births, deaths, merges, and splits of features; visualization algorithms for tracking graphs

## Dependencies

FTK requires CMake to build the libraries and executables.  Optional dependencies include ParaView (>=5.8.0 recommended), VTK, Qt5, MPI, netCDF, parallel-netcdf, HDF5, ADIOS2, MPSolve, and CUDA.

## Usage

You may use FTK as ParaView plugins, Python bindings, or command line interface.   

### FTK ParaView Plugins

FTK provides ParaView plugins to allow users track critical points (maxima, minima, and saddles) in scalar field data.  In order to build the plugins, we recommend to build and use 

```bash
$ git clone https://github.com/hguo/ftk $FTK_SOURCE_DIR
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_PARAVIEW=ON -DParaView_DIR=$YOUR_ParaView_Build
$ make
```

If built successfully, you will see the plugins binary as `lib/paraview-5.8/plugins/FTK/FTK.so`.  Open the "Plugin Manager" in ParaView, and load this binary with "Load New..." button, and then select and load FTK in the list.  

We demo the use the 2D critical point tracking filter (`vtkCriticalPoint2DTracker`) with a dataset.  The input of this filter must be a 3D volumetric data that stacks 2D time-varying scalar fields in the Z direction.  In this demo, we first add a synthetic 3D volume data by using Sources / FTK / Spiral2DSource.  We then track the trajectories of 2D critical points with Filters / FTK / CriticalPoint2DTracker.  The output trajectires can be visualized as tubes and color-coded by their types, scalar values, or IDs.  In this demo, the time-varying scalar field is defined in closed form: 

$f(x,y,t)=cos(x\cos t - y\sin t) \sin(x\sin t + y\cos t),$

where $x$ and $y$ are 2D coordinates and $t$ is time.  We discretize the $x,y$ domain into a $128\times 128$ regular grid and the time domain into 10 timesteps.  Local maximum are defined as the loci of points that $(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial x})=0$ and both eigenvalues of the Hessian of $f$ (in terms of $x$ and $y$) are negative.  We use a sweep-and-trace algorithm to first localize local maximum and trace the maximum over space-time.  We first mesh the scalar field with a 3D regular simplex mesh and check every 2-elements (faces) meets the criteria.  We then do the connected component labeling; two faces are connected if each of them has a local maxima and share the same 3-element (tetrahedra).  The trajectories are then constructured from the connected components.  

### FTK Python Bindings (PyFTK)

FTK Python bindings is based on [pybind11](https://github.com/pybind/pybind11).  

#### Build PyFTK

##### Build PyFTK with `setuptools`

Notice that CMake is required to build PyFTK.  Advanced build options is currently not possible to configure with `setuptools`.  Refer to the next section to build PyFTK directly with CMake.

```bash
$ cd $FTK_SOURCE_DIR
$ python setup.py install
```

##### Build PyFTK with CMake

```bash
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_PYFTK=ON
$ make
```

The output PyFTK binary will be in the `lib` directory.

#### PyFTK Usage

PyFTK provides synthetic data generators and critical point trackers.  The following is an example of tracking critical points in a synthetic spiral woven data:

```python
import pyftk

data = pyftk.synth.spiral_woven(10, 10, 20) # generate synthetic spiral woven data (width=10, height=10, and 20 timesteps).  The shape of data is (10, 10, 20)
result = pyftk.track_critical_points_2d_scalar(data) # track critical points in a scalar field
print(result)
```

The results are trajectories organized in a list: 

```
[{'length': 9, 'trace': [{'x': 2.275077079338536, 'y': 2.0, 't': 2.843946435964648, 'type': 'min', 'scalar': -0.7349697808320285}, {'x': 2.3009922790096073, 'y': 2.057205556154771, 't': 3.0, 'type': 'min', 'scalar': -0.7126261556354363}, {'x': 2.316376550504984, 'y': 2.0789601019629704, 't': 3.0789601019629704, 'type': 'min', 'scalar': -0.6994583185227987}, {'x': 2.3396684290296013, 'y': 2.109042720626548, 't': 3.339668429029601, 'type': 'min', 'scalar': -0.6203974444741183}, {'x': 2.4602960605411885, 'y': 2.367439624426215, 't': 4.0, 'type': 'min', 'scalar': -0.502426092806519}, {'x': 2.5836144734591056, 'y': 2.5204553926376145, 't': 4.520455392637614, 'type': 'saddle', 'scalar': -0.3968294787319291}, {'x': 2.587217124155211, 'y': 2.5205274563826645, 't': 4.587217124155211, 'type': 'saddle', 'scalar': -0.37723450315450113}, ...
```

### FTK C++ Libraries

You may use FTK as a C++ library.  The installation will also generate FTKConfig.cmake in the installation path, such that you can use `find_package(FTK)` to find and use FTK in your CMakeLists.txt

```bash
$ git clone https://github.com/hguo/ftk $FTK_SOURCE_DIR
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DCMAKE_INSTALL_PREFIX=$FTK_INSTALL_DIR
$ make install
```

The installed files are organized as follows: 

```bash
$ tree $FTK_INSTALL_DIR
.
├── include
│   ├── ftk
│   │   ├── algorithms
│   │   │   ├── bfs.hh
...
│   └── hypermesh
│       ├── ndarray.hh
│       ├── regular_mesh.hh
...
└── lib
    └── cmake
        └── FTKConfig.cmake
```

#### Include FTK in your CMake project

You may use the FTK installation in your own CMakeLists.txt file:

```cmake
find_package(FTK REQUIRED)
include_directories (${FTK_INCLUDE_DIR})
```

When you configure your build, please specify FTK_DIR with CMake: 

```bash
$ cmake -DFTK_DIR=$FTK_INSTALL_DIR/lib/cmake
```

### Standalone executables

FTK currently provides three example executables, including 2D critical point tracking, 3D critical point tracking, and a demo for the tracking graph.

#### Build FTK examples

To build FTK examples without any dependences, please follow the following configuration:

```bash
$ cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_EXAMPLES=1
$ make
```

To build with Qt5, VTK, and NetCDF, you need additional arguments for cmake:

```bash
$ cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_EXAMPLES=1 -DCMAKE_PREFIX_PATH="$path_to_qt5/lib/cmake;$path_to_vtk/lib/cmake" -DNETCDF_DIR="$path_to_netcdf"
$ make
```

#### Run FTK examples

The followings are the guide to run the 2D critical point tracking example.  We use the same synthetic data in the ParaView instructions, run the example with following command, and the executable prints the trajectory:

```bash
$ ./examples/critical_point_tracking_2d/ex_critical_point_tracking_2d
We found 16 trajectories:
--Curve 0:
---x=(2.000000, 26.972244), t=5.869000, val=0.955384
---x=(2.056124, 26.857716), t=5.857716, val=0.958166
...
--Curve 15:
---x=(125.000000, 122.414673), t=1.241837, val=0.884089
---x=(124.709686, 122.709686), t=1.187109, val=0.912740
...
```

If you built the example with Qt, you can use the Qt/OpenGL based window to explore the trajectory.  You can rotate the space with mouse and increase/decrease timesteps with the left/right keys.  The color encodes the ID of different trajectories.  Please use `—qt` argument to run the executable, and the example screenshot is as follows.

```bash
$ ./examples/critical_point_tracking_2d/ex_critical_point_tracking_2d --qt
```

![](./docs/images/critical_point_tracking_2d.png)

If you built the example with VTK, you may also explore the trajectories with a VTK window by adding `—vtk` argument: 

```bash
$ ./examples/critical_point_tracking_2d/ex_critical_point_tracking_2d --vtk
```

### Use MPI

You may use MPI and CUDA to accelerate feature tracking with both distributed- and GPU-parallelism.  To build FTK with MPI, you need to use MPI C/C++ compilers: 

```bash
$ CC=mpicc CXX=mpicxx cmake ..
```

In order to build FTK with CUDA, you also need to specify the path to the CUDA installation:

```bash
$ cmake .. -DCUDA_TOOLKIT_ROOT_DIR=$YOUR_CUDA_TOOLKIT_DIR
```

To run the above example with both MPI and CUDA, use `mpiexec` and `-x cuda`

```bash
$ mpiexec -n $NUM_PROCS ./examples/critical_point_tracking_2d/ex_critical_point_tracking_2d -x cuda
```

## Applications that use FTK

* [vortexfinder2](https://github.com/hguo/vortexfinder2): Vortex finder for time-dependent Ginzburg-Landau superconductor simulation data
