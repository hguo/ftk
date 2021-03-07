# FTK command line interface

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