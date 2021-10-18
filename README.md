This is the NUFEB development repository

NUFEB is an open source tool for modelling and simulating microbial communities.
The tool is based on the Individual-based Modelling (IbM) approach, 
where microbes are represented as discrete units and their
behaviour changes over time due to a variety of processes. 

NUFEB is built on top of the molecular dynamics simulator
LAMMPS (as a LAMMPS' user package), extended with IbM features. 
A wide range of biological, physical and
chemical processes are implemented to explicitly model microbial systems. 

NUFEB is a freely-available open-source code, distributed under the terms
of the GNU Public License.

NUFEB development has been funded by the UK’s EPSRC EP/K039083/1 
Newcastle University Frontiers in Engineering Biology (NUFEB) project.

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
lammps                  LAMMPS source code
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
git clone https://github.com/nufeb/NUFEB --recursive
</pre>

Build NUFEB:
<pre>
./install.sh
</pre>

### Running
Run a case in /examples after building NUFEB, for example:
<pre>
cd examples/biofilm-het
mpirun -np 4 ../../lammps/src/lmp_nufeb -in Inputscript.lammps
</pre>

---------------------------------------------------------------------------
Developers:

Bowen Li: bowen.li2@newcastle.ac.uk

Denis Taniguchi: denis.taniguchi@newcastle.ac.uk

