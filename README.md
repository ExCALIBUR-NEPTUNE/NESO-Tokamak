# NESO-Tokamak
A NESO-based solver for modelling MAST-U.

## Building
The easiest way to build the app is to use the [Spack package manager](https://spack.readthedocs.io/en/latest/index.html).
The build process has been tested with **Spack v0.24.0**; later versions may also work, but aren't recommended.

1. Install Spack v0.24.0, via the [official instructions](https://spack.readthedocs.io/en/latest/getting_started.html#installation).
On Debian-based systems (e.g. Ubuntu) the following should work:

```bash
$ git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git

# Optionally modify .bashrc such that spack is always initialised in new shells
echo 'export SPACK_ROOT=$HOME/.spack' >> $HOME/.bashrc
echo 'source $SPACK_ROOT/share/spack/setup-env.sh' >> $HOME/.bashrc
export SPACK_ROOT=$HOME/.spack
source $SPACK_ROOT/share/spack/setup-env.sh
```
2. Install intel compilers
```bash
spack install intel-oneapi-compilers
```
3. Build NESO
```bash
git clone https://github.com/ExCALIBUR-NEPTUNE/NESO.git
cd NESO
git submodule update --init
. activate
spack install
deactivate
```

4. Build and clone NESO-Tokamak:
```bash
cd ..
git clone https://github.com/ExCALIBUR-NEPTUNE/NESO-Tokamak.git
cd NESO-Tokamak
spack env activate ./spack -p
spack install
```

Note that some of the dependencies (particularly nektar++) can take some time to install.

## Running examples
Configuration files for 2D and 3D examples can be found in the `./examples` directory.
An easy way to run the examples is (from the top level of the repository):
```bash
# Load the mpich spack module
spack load mpich
# Set up and run an example via the helper script
./scripts/run_eg.sh [example_name] [example_dimension]
```

- This will copy `./examples/[example_name]/[example_dimension]` to `./runs/[example_name]/[example_dimension]` and run the solver in that directory with 4 MPI processes by default.
- To run with a different number of MPI processes, use `<-n num_MPI> `
- The solver executable is assumed to be in the most recently modified `spack-build*` directory. To choose a different build location, use `<-b build_dir_path>`.

To change the number of openMP threads used by each process, use
```bash
# Run an example with OMP and MPI
OMP_NUM_THREADS=[nthreads] ./scripts/run_eg.sh [example_name]
```

## Postprocessing
- Particle output is generated in an hdf5 file. See [particle_velhist.py](./postproc/particle_velhist.py) for an example of reading and postprocessing (requires [these packages](./postproc/requirements.txt)).
- Fluid output is generated as nektar++ checkpoint files. One way to visualise is them is to convert to .vtu using Nektar's `FieldConvert` tool and then use Paraview:
```bash
# Convert to vtu using FieldConvert (requires xml files as args)
cd runs/[example_name]/[example_dimension]
spack load nektar
FieldConvert [config_xml] [mesh_xml] [chk_name] [vtu_name]
# e.g. FieldConvert single_field.xml mastu.xml single_field_100.chk single_field.vtu
``` 
