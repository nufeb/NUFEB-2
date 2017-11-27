This is NUFEF development repository.

- doc: NUFEB user manuals and technical report for current and old NUFEB versions

- gui: a graphical user interface for NUFEB simulation (not compatible with the latest version)

- examples: IBm running cases 

- lammps5Nov16: LAMMPS code and NUFEB package

- openfoam-coupling: code for coupling LAMMPS with OpenFoam

- post-process: code for visualizing simulation result

- release: released version, see the public repository https://github.com/nufeb

- unittests: unit tests for NUFEB-1.0 functions

### How to build and run NUFEB simulation

This quick guide provides information about how to build and run the code

NUFEB requires GCC/G++ for a successful build. 

To build this code, first move to the /lammps5Nov16/src/STUBS directory, and make a dummy MPI lib:

$ make

Then, install the NUFEB and granular packages in /lammps5Nov16/src directory with the following commands:

$ make yes-USER-NUFEB

$ make yes-GRANULAR

In /lammps5Nov16/src directory, execute the following
 command to build the NUFEB executable:

$ make serial

Now, you can run the simulation by going to one of the subdirectoies in /inputs/examples/ directory and run:

$  ../../lammps5Nov16/src/lmp_serial < Inputscript.lammps
