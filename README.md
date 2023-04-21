<p align="center">
  <img src="https://github.com/nufeb/NUFEB-dev/blob/master/docs/source/images/log.png" | width=500/>
</p>

<h1 align="center">NUFEB: Agent-based Simulator for Microbial System Modelling</h1>


This is the NUFEB development repository.

NUFEB (Newcastle University Frontier in Engineering Biology) is an open source tool for 3D agent-based simulation of microbial systems.
The tool is built on top of the molecular dynamic simulator [LAMMPS](https://lammps.sandia.gov), and extended with features for microbial modelling. 
For more details about NUFEB and its agent-based model, please refer to our [software paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007125).

NUFEB is distributed under the terms of the GNU Public License. The development has been funded by the UK’s EPSRC EP/K039083/1 Frontiers in Engineering Biology project.

Online user manual is avaiable at: https://nufeb.readthedocs.io/en/master/ 

<p align="center">
  <img src="https://github.com/nufeb/NUFEB-dev/blob/master/docs/source/images/IBM-NUFEB.png">
</p>

---------------------------------------------------------------------------

The NUFEB distribution includes the following files and directories:
<pre>
README                  this file 
LICENSE                 the GNU General Public License (GPL)
install.sh              script for building NUFEB 
uninstall.sh            script for uninstalling NUFEB 
doc                     user manual and other documentation 
examples                test problems and cases used in publications 
lib                     libraries NUFEB can be linked with 
lammps_stable_29Oct2020 LAMMPS source code
src                     source files 
thirdparty              thirdparty tools
</pre>

NUFEB requires cmake, git, gcc/g++, openmpi, png-dev and ffmpeg for a successful build.
You can run the following commands to install the packages, 

On Ubuntu:
<pre>
sudo apt update
sudo apt-get install cmake git-core g++ openmpi-bin openmpi-common libopenmpi-dev libpng-dev ffmpeg
</pre>

On CentOS:
<pre>
sudo yum update
sudo yum install cmake git gcc-c++ openmpi openmpi-devel libpng-dev ffmpeg
</pre>

### Getting source
Get NUFEB source code and submodules for its thirdparty libraries:
<pre>
git clone https://github.com/nufeb/NUFEB-dev
</pre>

Build NUFEB:
<pre>
./install.sh
</pre>

### Running
Run a case in /examples after building NUFEB, for example:
<pre>
cd examples/biofilm-heterotroph
mpirun -np 4 ../../nufeb_mpi -in Inputscript.lammps
</pre>

---------------------------------------------------------------------------
Developers:

Bowen Li: bowen.li2@newcastle.ac.uk

Denis Taniguchi: denis.taniguchi@newcastle.ac.uk

Joseph E. Weaver: joe.weaver@newcastle.ac.uk
