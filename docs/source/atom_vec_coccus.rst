.. index:: atom_style coccus

atom_style coccus command
==========================

Syntax
""""""

.. parsed-literal::

    atom_style coccus

Description
""""""""""""""

Define the *coccus* style of atoms in a simulation. 
This command must be used before a simulation is setup via a
`read_data* <https://docs.lammps.org/read_data.html>`_, 
`read_restart* <https://docs.lammps.org/read_restart.html>`_, or
`create_box* <https://docs.lammps.org/create_box.html>`_ command.

.. image:: images/coccus.png
   :scale: 25% 
   :align: center
   
For the *coccus* style, atoms are represented as spheres that model
spherical-shaped microbes.
Each microbe stores a set per-atom attributes,
including *mass*, *biomass*, *outer_mass*, *diameter*, *outer_diameter*, *coordinate*,
as well as mechanical attributes used in physical (DEM) processes (e.g, *force*, *velocity*, etc).
 
The *diameter* and *outer_diameter* attributes specify the inner and outer sizes of a spherical microbe.
Some microorganism excrete extracellular polymeric
substances (EPS) which is initially accumulated as an extra shell beyond the microbes.
The *outer_diameter* is defined as the sum of EPS shell depth and the inner diameter,
while the *outer_mass* represents the weight of EPS shell.
The *mass* and *biomass* attributes correspond to the wet and dry weights of microbes, respectively.

Initial microbes and their attributes can be specified in 3 ways:

* use `read_data* <https://docs.lammps.org/read_data.html>`_ command to explicitly create each individual microbe with the initial attributes from a data file;
* use `create_atom* <https://docs.lammps.org/create_atom.html>`_ command to create microbes on a lattice, or a single microbe, or a random collection of microbes;
* use `read_restart* <https://docs.lammps.org/read_restart.html>`_ command to read previously saved system configuration from a restart file.
