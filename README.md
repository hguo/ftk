# FTK: The Feature Tracking Kit

[![Build Status](https://travis-ci.org/hguo/ftk.svg?branch=master)](https://travis-ci.org/hguo/ftk)

FTK is a library that scales, simplifies, and delivers feature tracking algorithms for scientific datasets.  You may use FTK as ParaView plugins, Python bindings, or command line interface.   See [arXiv:2011.08697](https://arxiv.org/abs/2011.08697) for more details on our feature tracking algorithms.

![](docs/images/critical_point_tracking_2d_paraview.png)

## Installation

To install FTK, use `spack install ftk`, `pip install pyftk`, or build from source; see [this page](docs/install.md) for more details on dependencies, compilation, and installation.  

## FTK command line interface

FTK provides one single executable `ftk`.  Follow the help information to explore all options the FTK CLI.

```bash
$ ftk --help
Usage:
  ./bin/ftk [OPTION...]
  -f, --feature arg             Feature type: critical_point, isosurface,
                                tdgl_vortex, or connected_component)
  -i, --input arg               Input file name pattern: a single file or a
                                series of file, e.g. 'scalar.raw',
                                'cm1out_000*.nc'
      --input-format arg        Input file format
...
```

#### Specify feature type (`-f`)

It is mandatory to use the argument `-f` specifies the type of features that are tracked, e.g. `critical_point` (`cp` for short), `isosurface` (`iso` for short), and `tdgl_vortex` (`tdgl` for short).  


#### Specify file inputs (`-i`)

Use `--input` or `-i` to specify input files.  The input argument typically needs to contain wildcards in order to read a list of files of the time-varying data.  For example, assuming the working directory has a series of files, each of which contains one timestep of the 2D regular-grid scalar field data: `scalar-000.raw`, `scalar-001.raw`, ... `scalar-255.raw`; to track critical points: 

```bash
$ ftk -f cp --input "scalar-*.raw" --width 128 --height 128 --input-format float64 --output trajs.vtp
```

Note that it is highly recommended to **use quotes** for the `-i` argument; otherwise the executable may only read the first file, and the behavior may very depending on which shell you use.  It is also important to **make sure that wildcards lead to the correct sequence** of files.  A counterexample is `0.raw`, `1.raw`, ... `10.raw`, `11.raw`, ..., `100.raw`; in this case, the sequence could be parsed as 0, 1, 10, 100, 11, ..., 2, 20, 21, which is undesired and could lead to meaningless results.  


##### Input format

Use `--input-format` to explicitly specify the input format, if the format is not possible to be determined based on file extensions.  Acceptable arguments for `--input-format` include `float32`, `float64`, `vti` (VTK image data format), `h5` (HDF5), and `nc` (netCDF).  In case file names end with `.vti`, `.nc`, and `.h5`, one may omit `--input-format`.  

The data of `float32` and `float64` are in the "block-of-value" format, and we assume the data are stored in the row-major order, a,k,a the C-style.  If such data contain mmultiple components (e.g. u and v), the dimension of components are the most fast-changing dimension in the file, e.g. (u0, v0, u1, v1, ...).  One needs to reorder the data if they are not in the correct order.  See the "input variables" section for more information.


##### Input dimensions

It is mandatory to specify `--width`, `--height`, `[--depth]` if the inputs are in `float32` or `float64` formats.  Note that `--depth` only applies to 3D data.

##### Input timesteps

We assume the number of timesteps are the number of files (unless in the case that each netCDF file contains multiple timesteps).  One may override the number of timesteps by explicitly use the `--timesteps` option, but the number must be less than or equals to the number of available timesteps. 

##### Input variables

Depending on file formats and number of components, it may or may not be necessary to specify `--var`.  If one single variable (e.g. `scalar`) is used, specify the variable name (e.g. `--var scalar`); if multiple variables are used, use comma between variable names, e.g. `--var u,v,w` for `u`, `v`, and `w`. 

- `float32`/`float64` (single component): not necessary
- `float32`/`float64` (multi component): necessary, e.g. `--var u,v`; in this case, components, width, height, depths are presumed to be stored in the raw-major order
- `nc` and `h5`: necessary regardless of single- or multi-components; the dimension of each variable must be identical
- `vti`: recommended but not necessary; the default variable will be used if `--var` is not specified


#### Use synthetic inputs in lieu of file inputs

See [this page](docs/synthetic.md) for more details on synthetic inputs.  Use `--synthetic` to use synthetic data instead of file inputs for demonstration and testing purposes.  Available options include `woven`, `double_gyre_2d`, `merger_2d`, `moving_extremum_2d`, `moving_extremum_3d`.  For example, to track critical in the 2D woven synthetic data:

```bash
$ ftk -f cp --synthetic woven --output woven.txt 
```

The dimensions and number of timesteps, which default values for each synthetic case, may be overridden by `--width`, `--height`, `--depth`, and `--timesteps` options.  Note that the `--depth` option is not applicable for 2D synthetic cases.  For example, to track critical points in a 30x30x30x10 3D moving extremum case:

```bash
$ ftk -f cp --synthetic momving_extremum_3d --width 30 --height 30 --depth 30 --timesteps 10 --output woven.txt 
```

#### Outputs

Use `--output` to specify output file name(s).  Outputs are usually in one single file; in cases that outputs are in multiple files, use wildcards.  For example, `out-%03d.vtp` will lead to a sequence of files (`out-000.vtp`, `out-001.vtp`, ...).  

##### Output formats

Possible output formats are plain text (`.txt`), JSON (`.json`), `vtkPolyData` (`.vtp`), `vtkUnstructuredGrid` (`.vtu`), and binary (for debugging purposes only), depending on the feature type and whether FTK is compiled with external dependencies.  The executable will automatically determine the output format, unless `--output-format` is specified. 

##### Output types

In general, there are four types of outputs, specified by the `--output-type` option:

- `traced` (default).  Outputs are trajectories of features in one single file, e.g. critical point trajectories (in text, JSON, or `.vtp` formats), surface trajectories of vortex lines (in `.vtp` or `.vtu` formats), and volume trajectories of isosurfaces (in `.vtu` format)
- `sliced`.  Outputs are features in individual timesteps, each timestep corresponds to an output file, e.g., critical points (in text, JSON, or `.vtp` formats), vortex lines (in `.vtp` format), and isosurfaces (in `.vtp` format)
- `discrete`.  The single output file contains untraced feature points in spacetime.  
- `intercepted`.  A series of output files; each is a subset of traced features for the given duration of time (specified by `--intercept-length`), currently only available to critical point tracking.


#### Parallel execution

Use `mpiexec` for distributed execution; use `--nthreads` and `--accelerator` to specify number of threads and if a hardware accelerator (GPU) is used. 

##### MPI

The FTK executable recognizes the option to use multiple processes if `mpiexec` is used.  It is recommended to carefully specify `--nthreads` to avoid over-subscribing resources if multiple processes are on the same computing node:

```bash
$ mpiexec -n $NUM_PROCS ftk --nthreads $NUM_THREADS_PER_PROC ...
```

##### POSIX Threads

By default, the executable uses the maximum number of hardware threads, unless `--nthreads` is specified.  

##### CUDA

Use `--accelerator cuda` if FTK is compiled with CUDA and an NVIDIA GPU is available. 

## FTK for ParaView

### Building ParaView plugins

FTK provides ParaView plugins to allow users track critical points (maxima, minima, and saddles) in scalar field data.  In order to build the plugins, we recommend to build and use 

```bash
$ git clone https://github.com/hguo/ftk $FTK_SOURCE_DIR
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DFTK_BUILD_PARAVIEW=ON -DParaView_DIR=$YOUR_ParaView_Build
$ make
```

If built successfully, you will see the plugins binary as `lib/paraview-5.8/plugins/FTK/FTK.so`.  Open the "Plugin Manager" in ParaView, and load this binary with "Load New..." button, and then select and load FTK in the list.  To check if ParaView plugins are correctly built by reproducing the results in the above figure, use "Sources-->FTK-->SpiralWoven2DSource", "Filters-->FTK-->CriticalPointTracker2D",followed by the "Tube" filter in ParaView.

### Using ParaView plugins

We demonstrate the use the 2D critical point tracking filter (`vtkCriticalPoint2DTracker`) with a dataset.  The input of this filter must be a 3D volumetric data that stacks 2D time-varying scalar fields in the Z direction.  In this demo, we first add a synthetic 3D volume data by using Sources / FTK / Spiral2DSource.  We then track the trajectories of 2D critical points with Filters / FTK / CriticalPoint2DTracker.  The output trajectires can be visualized as tubes and color-coded by their types, scalar values, or IDs.  In this demo, the time-varying scalar field is defined in closed form: 

$f(x,y,t)=cos(x\cos t - y\sin t) \sin(x\sin t + y\cos t),$

where $x$ and $y$ are 2D coordinates and $t$ is time.  We discretize the $x,y$ domain into a $128\times 128$ regular grid and the time domain into 10 timesteps.  Local maximum are defined as the loci of points that $(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial x})=0$ and both eigenvalues of the Hessian of $f$ (in terms of $x$ and $y$) are negative.  We use a sweep-and-trace algorithm to first localize local maximum and trace the maximum over space-time.  We first mesh the scalar field with a 3D regular simplex mesh and check every 2-elements (faces) meets the criteria.  We then do the connected component labeling; two faces are connected if each of them has a local maxima and share the same 3-element (tetrahedra).  The trajectories are then constructured from the connected components.  

## FTK for Python (PyFTK)

You can install PyFTK with `pip`.  The only dependency in the current release is `numpy`.  If you would like to build PyFTK from source and use advanced build options, see [this page](docs/install.md).

```bash
$ pip3 install pyftk
```

PyFTK provides synthetic data generators (`pyftk.synthesizers`), feature extractors (`pyftk.extractors`), and feature trackers (`pyftk.trackers`).  Currently, PyFTK only supports critical points.  The following is an example of tracking critical points in a synthetic spiral woven data:

```python
>>> import pyftk
>>> data = pyftk.synthesizers.spiral_woven(10, 10, 20) # generate synthetic spiral woven data (width=10, height=10, and 20 timesteps).  The shape of data is (1, 10, 10, 20)
>>> result = pyftk.trackers.track_critical_points_2d_scalar(data) # track critical points in a scalar field
>>> print(result)
```

The results are trajectories organized in a list: 

```
[{'length': 9, 'trace': [{'x': 2.275077079338536, 'y': 2.0, 't': 2.843946435964648, 'type': 'min', 'scalar': -0.7349697808320285}, {'x': 2.3009922790096073, 'y': 2.057205556154771, 't': 3.0, 'type': 'min', 'scalar': -0.7126261556354363}, {'x': 2.316376550504984, 'y': 2.0789601019629704, 't': 3.0789601019629704, 'type': 'min', 'scalar': -0.6994583185227987}, {'x': 2.3396684290296013, 'y': 2.109042720626548, 't': 3.339668429029601, 'type': 'min', 'scalar': -0.6203974444741183}, ...
```


## FTK C++ Libraries

You may use FTK as a C++ library.  The installation will also generate FTKConfig.cmake in the installation path, such that you can use `find_package(FTK)` to find and use FTK in your CMakeLists.txt

```bash
$ git clone https://github.com/hguo/ftk $FTK_SOURCE_DIR
$ mkdir $FTK_SOURCE_DIR/build && cd $FTK_SOURCE_DIR/build
$ cmake .. -DCMAKE_INSTALL_PREFIX=$FTK_INSTALL_DIR
$ make install
```

#### Including FTK in your CMake project

You may use the FTK installation in your own CMakeLists.txt file:

```cmake
find_package(FTK REQUIRED)
include_directories (${FTK_INCLUDE_DIR})
target_link_library (${YOUR_TARGET} FTK::libftk)
```

Specify FTK_DIR with CMake: 

```bash
$ cmake -DFTK_DIR=$FTK_INSTALL_DIR/lib/cmake
```

### FTK library components

* Hypermesh: data structures for high-dimensional meshes and mesh elements including *n*-simplices, *n*-cubes, and *n*-prisms; utilities to generalize given 2D/3D structured/unstructured meshes into 3D/4D spacetime meshes

![](./docs/images/regular_simplex_subdivision.svg)

* Numeric: root-find algorithms for inverse interpolations and parallel vector operators in *n*-simplices, *n*-cubes, and simplex-prisms; lightweight linear algebra utilities to support root-finding
* CCL: connected component labeling algorithm for building feature tracking algorithms
* Geometry: utilities to transform connect components to geometry for visualization and analysis
* Tracking graph: data structures to record births, deaths, merges, and splits of features; visualization algorithms for tracking graphs

## Contact

* [Hanqi Guo](http://www.mcs.anl.gov/~hguo/), [hguo@anl.gov](mailto:hguo@anl.gov)

## Publications

* Hanqi Guo, David Lenz, Jiayi Xu, Xin Liang, Wenbin He, Iulian R. Grindeanu, Han-Wei Shen, Tom Peterka, Todd Munson, and Ian Foster, "FTK: A High-Dimensional Simplicial Meshing Framework for Robust and Scalable Feature Tracking" [arXiv:2011.08697](https://arxiv.org/abs/2011.08697) [cs.GR], 2020.

* Jiayi Xu, Hanqi Guo, Han-Wei Shen, Mukund Raj, Xueqiao Xu, Xueyun Wang, Zhehui Wang, and Tom Peterka, "Asynchronous and Load-Balanced Union-Find for Distributed and Parallel Scientific Data Visualization and Analysis." IEEE Transactions on Visualization and Computer Graphics (Proc. IEEE PacificVis 2021), [arXiv:2003.02351](https://arxiv.org/abs/2003.02351) [cs.DC], 2021. (Accepted)
