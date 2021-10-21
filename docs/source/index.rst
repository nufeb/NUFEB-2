.. NUFEB-manual documentation master file, created by
   sphinx-quickstart on Tue Oct 19 16:27:51 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NUFEB Documentation
========================================

NUFEB (Newcastle University Frontier in Engineering Biology) 
is an open source tool for 3D individual-based simulation of microbial communities. 

NUFEB is built on top of the molecular dynamic simulator `LAMMPS <https://lammps.sandia.gov/>`_ (version: stable_29Oct2020), 
and extended
with features for microbial modelling. 
As a result, NUFEB has the realistic explicit models of biological, chemical and
physical processes, including individual microbes in different shapes. 
The current implementation takes full advantage of the LAMMPS parallel infrastructure,
allowing for massive parallel simulations on both CPU- and GPU-based systems.

NUFEB is distributed under the terms of the GNU Public License. 
The development has been funded by the UK’s EPSRC EP/K039083/1 Frontiers in Engineering Biology project.

| Development team:
| Bowen Li, bowen.li2@newcastle.ac.uk
| Denis Taniguchi: denis.taniguchi@newcastle.ac.uk
| Joe Weaver: joe.weaver@newcastle.ac.uk

----------

This site contains user and programmer documentation for the NUFEB use. If you find
errors or omissions in the manual or have suggestions for useful information to add, please
send an email to the developers.

----------

.. _user_documentation:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: User Guide
   :name: userdoc
   :includehidden:
   
   intro
   install
   commands
   packages
   

.. _programmer_documentation:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: Programmer Guide
   :name: progdoc
   :includehidden:

   dev_new_growth


.. _command_reference:
.. toctree::
   :name: reference
   :maxdepth: 1
   :caption: Command Reference

   list_general
   list_microbe
   list_biology
   list_physics
   list_chemsitry
   list_compute
   list_post
   list_lammps
   
   
******************
Indices and tables
******************

.. only:: html

   * :ref:`genindex`
   * :ref:`search`