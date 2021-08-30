## Replicability with CLI

This document serves as a guideline to reproduce the result in Figure 14 of this [preprint](https://arxiv.org/abs/2011.08697) with FTK's CLI:

```bash
$ ftk --feature critical_point --synthetic merger_2d --output merger_2d.vtp
```

The above command tracks critical poitns in the synthetic merger data and writes critical point trajectories in the `merger_2d.vtp` file.  The `.vtp` output can be openned and visualized with [ParaView](https://www.paraview.org/).  

Note: 

* This example requires building FTK with VTK
* This example requires prebuilt ParaView to generate visualization results
* The FTK CLI does not produce any visualization
* This example does **NOT** require FTK built with ParaView

The following documents building FTK with VTK in the vanila macOS Mojave 10.14.6 installation with Apple clang version 11.0.0.  Assuming CMake is already installed (from source, Homebrew, or other package management tools).


### Building VTK

Download the VTK 9.0.1 source ([](https://vtk.org/download/)), and then build and install the VTK library (this may take a while):

```bash
$ tar zxf VTK-9.0.1.tar.gz
$ cd VTK-9.0.1
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/your/vtk/install/path
$ make && make install
```

### Building FTK with VTK

Assuming FTK source is already in `$FTK_SOURCE_DIR`:

```bash
$ cd $FTK_SOURCE_DIR
$ mkdir build && cd build
$ cmake .. -DFTK_USE_VTK=ON -DCMAKE_PREFIX_PATH=/your/vtk/install/path -DCMAKE_INSTALL_PREFIX=/your/ftk/install/path
$ make && make install
```

### Run the executable

After FTK is installed, add FTK and VTK binary/library paths to environmental variables before running the executable:

```bash
$ export PATH=/your/ftk/install/path:$PATH
$ export DYLD_LIBRARY_PATH=/your/vtk/install/path/lib:/your/ftk/install/path/lib:$DYLD_LIBRARY_PATH
```

Now, use the command in the begging of this document to generate `merger_2d.vtp`.
