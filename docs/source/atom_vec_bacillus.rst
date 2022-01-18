.. index:: atom_style bacillus

atom_style bacillus command
============================

Syntax
""""""

.. parsed-literal::

    atom_style bacillus

Description
""""""""""""""

Define the *bacillus* style of atoms in a simulation. 
This command must be used before a simulation is setup via a 
`read_data* <https://docs.lammps.org/read_data.html>`_, 
`read_restart* <https://docs.lammps.org/read_restart.html>`_, or
`create_box* <https://docs.lammps.org/create_box.html>`_ command.

.. image:: images/rod.png
   :scale: 30% 
   :align: center
   
For *bacillus* style, atoms are represented as rods which model
rod-shaped microbes. 
Each microbe stores a set per-atom attributes, 
including *length*, *diameter*, *mass*, *biomass*, *coordinate*,
as well as the mechanical attributes for DEM simulation 
(e.g, *force*, *velocity*, *inertia*, *angular momentum*, etc).
 
In NUFEB, rod is modeled as cylinider with hemispherical caps. 
Therefore, cell *length* is the height of the cylinder and 
cell width is the *diameter* of the hemispherical caps.
The *mass* and *biomass* are the wet and dry weights of microbes, respectively. 

Initial microbes and their attributes can be specified in 3 ways: 
1) use `read_data* <https://docs.lammps.org/read_data.html>`_ command to 
explicitly create each individual microbe with the initial attributes from a data file;
2) use `create_atom* <https://docs.lammps.org/create_atom.html>`_ command
to create microbes on a lattice, or a single microbe, or a random collection of microbes;
or 3) use `read_restart* <https://docs.lammps.org/read_restart.html>`_ command to read 
previously saved system configuration from a restart file.
