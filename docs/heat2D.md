# Example: In situ critical point tracking in a 2D heat transfer simulation

This example demonstrates the use of FTK command line to track scalar critical points in a heat transfer simulation in situ via ADIOS2.  To run this example, FTK needs to be built with ADIOS2 and HDF5.  There is no need to modify the simulation code.

## Build the simulation

Checkout [this](https://github.com/hguo/adiosvm/tree/master/Tutorial/heat2d/cpp) 2D heat transfer simulation code from ADIOS2 examples.  Edit `make.settings` to set `ADIOS2_DIR` to the same ADIOS2 installation that FTK is built with.  Build the simulation code with `make`.

## Edit the configuration

Edit the `adios2.xml` and change the `Engine` of `SimulationOutput` to `InSituMPI`, `InSituAnalysis`, or `SST` for in situ analysis.  One may remain this file unchanged to produce the outputs while writing the output data to a `.bp` file.  

## Execute the simulation and FTK executable

```bash
$ mpiexec -n 12 ./heatSimulation sim.bp 4 3 5 10 200 1 : -n 2 /path/to/ftk --adios-config adios2.xml --adios-name SimulationOutput -f cp --input sim.bp --var T --output-type traced --output heat.vtp --stream
```

This command invokes the simulation and FTK with 12 and 2 processes, respectively.  The data resolution is 20 by 30.  Trajectories of critical points are written into the `heat.vtp` file at the exit.

## Visualize the outputs

Open `heat.vtp` with ParaView to visualize the trajectories of critical points.
