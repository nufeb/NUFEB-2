Introduction
============

This section gives a brief overview of the NUFEB software and model,
lists NUFEB key features, 
and acknowledges people who
have contributed to NUFEB software and model developments.

.. contents:: 
		:local:
		:depth: 1
   




.. _intro_1:

What is NUFEB
----------------

NUFEB is an open-source individual based model tool
for simulating the 3D dynamics of cell populations at the micro-scale.
The tool is built on top of the molecular dynamic simulator LAMMPS,
extended with features for biological system modelling.
NUFEB offers the realistic explicit models of biological, chemical and
physical processes, as well as individual cell in different shapes.

NUFEB takes full advantage of the LAMMPS parallel infrastructure,
allowing for massive parallel simulations on both CPU- and GPU-based parallel computers. 
On the other hand, it can also be built and run on laptop or desktop machine.


What is Individual-based Model
--------------------------------

Individual-based Model (IbM) is a computational model for 
simulating the actions and interactions of autonomous individuals,
which, in our case, are the cells.

IbMs precisely capture how the heterogeneity of 
individual organisms and local interactions influence emergent behaviours, 
such as bacterial colony morphology.

In IbM, microbes are represented as rigid particles, 
each with their own set of biological and physical properties. 
These properties are affected by internal or external processes, 
resulting in cell growth, decay, motility, etc.

NUFEB general features
------------------------

* open-source distribution
* runs on a single processor or in parallel
* CPU parallelism support for all features (MPI, OpenMP)
* GPU parallelism support for partial features (Kokkos)
* runs from an input script
* syntax for defining and using variables and formulas
* syntax for looping over runs and breaking out of loops
* run one or multiple simulations simultaneously (in parallel) from one script
* output to VTK or HDF5 data format
* output to PNG or JPEG image files
* easy to extend with new features and functionality
* build as library, invoke by other programme

NUFEB model features
------------------------

* spherical, rod, and ellipsoid cell representations
* flexible definition of cell distributions
* 3d computational domain
* continuous field for solute or gas distribution
* biological processes (growth, division, plasmid replication, etc)
* chemical processes (nutrient diffusion, pH, energy, etc)
* DEM-based physical processes (repulsive force, adhesion, etc)
* OpenFOAM integration for simulating fluid dynamics (CFD-DEM coupling)

Acknowledgments
------------------------

NUFEB development has been funded by the following sources:

* UK EPSRC EP/K039083/1 Frontiers in Engineering Biology project
* UK EPSRC IAA (Impact Acceleration Accounts) award


NUFEB development team
------------------------

NUFEB is developed at Newcastle University, UK. 

| The current NUFEB developers are:
| * `Bowen Li <https://www.ncl.ac.uk/computing/staff/profile/bowenli2.html>`_ (bowen.li2@newcastle.ac.uk)
| * `Denis Taniguchi <https://www.ncl.ac.uk/engineering/staff/profile/denistaniguchi.html>`_ (denis.taniguchi@newcastle.ac.uk)
| * `Joseph E. Weaver <https://joeweaver.github.io/>`_ (joe.weaver@newcastle.ac.uk)

| Past developers include:
| * Prashant Gupta 
| * `Curtis Madsen <https://sites.bu.edu/ckmadsen//>`_

| Special thanks to the following people who contributed to the NUFEB (Ib) model development:
| * `Jayathilake Pahala Gedara <https://www.oncology.ox.ac.uk/team/jayathilake-pahala-gedara>`_
| * `Valentina Gogulancea <https://www.ulster.ac.uk/staff/v-gogulancea>`_
| * Rebeca Gonzalez-Cabaleiro 
| * `Jonathan Sakkos <https://www.jonathanksakkos.com/>`_


Citing NUFEB
------------------------
The following paper describes NUFEB (v1.0)’s functionalities and implementation details:

Li B, Taniguchi D, Gedara JP, Gogulancea V, Gonzalez-Cabaleiro R, et al. (2019) 
NUFEB: A massively parallel simulator for individual-based modelling of microbial communities. 
PLOS Computational Biology 15(12): e1007125. https://doi.org/10.1371/journal.pcbi.1007125




