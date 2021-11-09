Installation on Linux & Mac
================================

This section describes how to install and run NUFEB on Linux and Mac OS. 

.. contents:: 
		:local:
		:depth: 1
   




.. _install_1:

Pre-compilation instructions
--------------------------------

Before compiling NUFEB, please make sure you have installed with the 
following packages depending upon the operating system used:

 *   cmake (https://cmake.org/)
 *   git (https://git-scm.com/)
 *   gcc/g++ (https://gcc.gnu.org/)
 *   openmpi (https://www.open-mpi.org/)
 *   libpng (http://www.libpng.org/)
 *   ffmpeg (https://www.ffmpeg.org/)
 

You can run the following commands to install the packages,

On **Ubuntu**:

 .. parsed-literal::

   sudo apt update
   sudo apt-get install cmake git-core g++ openmpi-bin openmpi-common libopenmpi-dev libpng-dev ffmpeg
   
On **CentOS**:

 .. parsed-literal::
   sudo yum update
   sudo yum install cmake git gcc-c++ openmpi openmpi-devel libpng-devel ffmpeg
   
On **MacOS**, you can use Homebrew to install the required libraries:

 .. parsed-literal::
   ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
 
then type the following command for the library installation:

 .. parsed-literal::
   brew install cmake git gcc open-mpi libpng ffmeg
   
   
Downloading NUFEB
--------------------------------

Use GIT checkout and update commands to get the NUFEB files once and then stay current. 
To do this, use the clone command to create a local copy of NUFEB Github repository:

 .. parsed-literal::
   git clone --recursive https://github.com/nufeb/NUFEB-dev.git
   
Once the command completes, a new directory named "NUFEB-dev" will be 
created on your machine which contains the latest NUFEB source code, examples, 
LAMMPS and thirdparty libraries. After initial cloning, 
as bug fixes and new features are added to NUFEB, 
you can stay up-to-date by typing the following GIT commands in the *NUFEB-dev* directory:

 .. parsed-literal::
   git checkout master
   git pull


Compiling NUFEB
--------------------------------

Use the script file ``install.sh`` provided in *NUFEB-dev* directory to compile the code:

 .. parsed-literal::
   ./install.sh [--options]
   
The [-\-options] 
is the additional settings of the compilation. The table below lists all available options:   

+--------------------+------------------------------------------------------------------------+
| **Options**        | **Description**                                                        |
+--------------------+------------------------------------------------------------------------+
| default            | mpi version (run on multiple CPUs)                                     |
+--------------------+------------------------------------------------------------------------+
| -\-serial          | serial version (run on single CPU only)                                |
+--------------------+------------------------------------------------------------------------+
| -\-gpu             | GPU version (support GPU parallelisation)                              |
+--------------------+------------------------------------------------------------------------+
| -\-enable-vtk  \*  | compile with VTK library to allow dumping vtk data files               |
+--------------------+------------------------------------------------------------------------+
| -\-enable-hdf5 \*  | compile with HDF5 library to allow dumping hdf5 data files             |
+--------------------+------------------------------------------------------------------------+
| -\-enable-plasmid  | compile with PLASMID optional package                                  |
+--------------------+------------------------------------------------------------------------+
| -\-static          | compile as a static library                                            |
+--------------------+------------------------------------------------------------------------+
| -\-shared          | compile as a shared library                                            |
+--------------------+------------------------------------------------------------------------+

\* 
Those options require external library located in *NUFEB/thirdparty* directory 
to be installed prior to the NUFEB compilation. 
You can install the library by running the corresponding script file, 
for example:

 .. parsed-literal::
	cd thirdparty
	./install-vtk.sh


It is possible to have more than one options. For example, running the command

 .. parsed-literal::
   ./install.sh --enable-vtk --enable-hdf5
   
will allow NUFEB simulation to output both vtk and hdf5 data formats.

When the installation finished, you should have an executable ``nufeb_mpi`` or
``nufeb_serial`` or ``nufeb_gpu`` in 
*NUFEB-dev* directory deponeding on configuration. The path to the executable 
will be automatically added to the system (.bashrc).

 .. note::
   For the sake of convenience, the executables built from install.sh are limited to mpi, gpu (cuda + mpi), and serial versions.
   More building options can be achieved by using traditional makefiles, see `LAMMPS user manual <https://docs.lammps.org/Build.html/>`_. 
   for the details.
   
   

Running NUFEB
--------------------------------

By default, NUFEB runs by reading commands from standard input. Thus if you run the NUFEB executable by itself, e.g.

 .. parsed-literal::
   nufeb_mpi

it will simply wait, expecting commands from the keyboard. 

Typically you should put commands in an input script and use I/O redirection, e.g.

 .. parsed-literal::
   mpirun -np 4 nufeb_mpi -in Inputscript.lammps
  
Here we are using parallel environments (MPI), 
using `mpirun -np X` command-line switch to specify the numbers of CPUs to run the simulation. 

Example cases are avaiable in *NUFEB-dev/examples* directory. To run one of them, 
go to its subdirectories and run NUFEB executable by passing in the input file, for example:

 .. parsed-literal::
  cd NUFEB-dev/examples/biofilm-het/  
  mpirun -np 4 nufeb_mpi -in Inputscript.lammps

or

 .. parsed-literal::
  cd NUFEB-dev/examples/biofilm-het/
  mpirun -np 4 nufeb_mpi -in Inputscript-vtk.lammps
  
if the VTK option is enabled. 
Output files will be generated and saved in a subdirectory during the simulation.
