.. index:: fix nufeb/property/cycletime

fix nufeb/property/cycletime command
=====================================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/property/cycletime 
    
* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to

Examples
""""""""

.. code-block:: 

   fix p_ctime HET nufeb/property/cycletime 
   
Description
"""""""""""

Assign cell cycle time attribute to each atom defined in *group-ID*.
The attribute tracks the elapsed time for its cell cycle.
Whenever a new atom is created or enters to the system
(e.g, through `create_atom* <https://docs.lammps.org/create_atom.html>`_
or :doc:`fix nufeb/divide/coccus <fix_divide_coccus>`),
the elapsed time for the cell cycle starts accumulating.
Once the cell completes the cell cycle phase through division,
the total elapsed time is store before reset to zero.

The pre-atom attribute can be dumped into various file formats (e.g, via `dump_vtk* <https://docs.lammps.org/dump_vtk.html>`_ command)