Run NUFEB
================================

This section describes how to run NUFEB 

.. contents:: 
		:local:
		:depth: 1
   

.. _run_nufeb:


Run NUFEB from terminal
--------------------------------

NUFEB run from the command line, reading commands from a file via the `-in` command line flag, or from standard input. 
NUFEB executable is located in the *NUFEB-dev* directory. 
The name of the executable is ``nufeb_mpi`` or ``nufeb_gpu`` or ``nufeb_serial`` depending on the selected option when compiling NUFEB.

You normally run the nufeb command in the directory where your input script is located. 
That is also where output files are produced by default, unless you provide specific other paths in your input script or on the command line.

NUFEB provides numbers of examples available in the *NUFEB-dev/examples* directory.
For example: 

 .. parsed-literal::
   cd NUFEB-dev/examples/biofilm-het
   mpirun -np 4 nufeb_mpi -in nufeb.in

The above commands will run the input script *nufeb.in* in the *NUFEB-dev/biofilm-het* directory using NUFEB MPI executable with 4 CPU processors.
`mpirun -np X` command-line switch to specify the numbers of processors to run the simulation.

 .. parsed-literal::
  cd NUFEB-dev/examples/biofilm-het/
  nufeb_serial -in nufeb.in

The above commands will run the *biofilm-het* example using NUFEB serial executable with single processor.

 .. parsed-literal::
  cd NUFEB-dev/examples/biofilm-het/  
  mpirun -np 2 nufeb_gpu -k on g 2 -sf kk -pk kokkos newton off neigh half binsize 2e-6 -in nufeb.in

The above commands will run the *biofilm-het* examaple using NUFEB GPU executable with 2 GPUs and 2 CPUs.
In particular, the `-k on` command-line switch enable Kokkos support; `g 2` specify the use of 2 GPUs to run the simulation.
`-sf kk` enables NUFEB to use Kokkos variant styles of any commands defined in the inputscript; 
`-pk kokkos newton off neigh half binsize 2e-6` re-define the setting in the inputscript.

You can find more information about `command line options* <https://docs.lammps.org/Run_options.html>`_ and
`screen and logfile output* <https://docs.lammps.org/Run_output.html>`_ from LAMMPS user manual.