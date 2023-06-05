.. index:: fix nufeb/property/generation

fix nufeb/property/generation command
=====================================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/property/generation
    
* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to

Examples
""""""""

.. code-block:: 

   fix p_gen HET nufeb/property/generation
   
Description
"""""""""""
Assign cell generation attribute to each atom defined in *group-ID*
The attribute tracks the cell generation count. The count is initialised to 1 whenever
a new cell is created or enters to the system
(e.g, through `create_atom* <https://docs.lammps.org/create_atom.html>`_
or :doc:`fix nufeb/divide/coccus <fix_divide_coccus>`).
Once the cell divide, each daughter cell has its generation count incremented.

The pre-atom attribute can be dumped into various file formats (e.g, via `dump_vtk* <https://docs.lammps.org/dump_vtk.html>`_ command)