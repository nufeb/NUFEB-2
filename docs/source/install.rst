Install NUFEB on Linux & Mac
================================

This section describes how to install NUFEB-2 on Linux and Mac OS.

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
   git clone https://github.com/nufeb/NUFEB-2.git
   
Once the command completes, a new directory named "NUFEB-2" will be 
created on your machine which contains the latest NUFEB source code, examples, 
LAMMPS and thirdparty libraries. After initial cloning, 
as bug fixes and new features are added to NUFEB, 
you can stay up-to-date by typing the following GIT commands in the ``NUFEB-2`` directory:

 .. parsed-literal::
   git checkout master
   git pull


Compiling NUFEB
--------------------------------

Use the script file ``install.sh`` provided in ``NUFEB-2`` directory to compile the code:

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
Those options require external libraries located in the ``NUFEB-2/thirdparty`` directory 
to be installed prior to NUFEB compilation. 
You can install the libraries by running the corresponding script file, 
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
``NUFEB-2`` directory deponeding on configuration. The path to the executable 
will be automatically added to the system (.bashrc).

 .. note::
   For convenience, the executables built from install.sh are limited to mpi, gpu (cuda + mpi), and serial versions.
   More building options can be achieved by using traditional makefiles, see `LAMMPS user manual <https://docs.lammps.org/Build.html/>`_
   for the details.
   
   

