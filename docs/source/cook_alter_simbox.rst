Altering the dimensions of the simulation box
=============================================


Problem
-------

You wish to adjust the width, length, or heigh of the simulation domain.

Solution
--------

Edit the appropriate dimensions on lines 6-8 of the ``atom.in`` file

.. code-block::
    :linenos:
    :emphasize-lines: 6,7,8

     NUFEB Simulation
    
         40 atoms
         2 atom types
  
     0.0e-04   1e-04  xlo xhi
     0.0e-04   1e-04  ylo yhi
     0.0e-04   5e-04  zlo zhi


Discussion
----------

NUFEB simulations often specify the intial microbial layout, including the simulation dimensions, in an external file by making use of the LAMMPS ``read_data`` command. This file is traditionally named ``atom.in``.  Lines 6-8 of this file specify the minimum and maximum bounds of the x, y, and z dimensions, with z being the vertical 'height' dimension.  The values can be given as either plan digits or in scientific notation, and the default SI unit is metres.

.. tip::
    Units are specified in the inpufile using the LAMMPS ``unit`` command
    It can be easy to mistakenly alter the wrong dimension, it is highly suggested to confirm them in a test run by visualizing results with paraview.
    
.. warning::

   Altering the simulation box might result in a situation where the diffusion grid voxels specified by the ``grid_style`` command in the inputscript may not evenly fit. This will cause an error at runtime and can be fixed by specifing a voxel size which evenly fits within the new dimensions.

.. seealso::
    Relevant section of the inputscript documentation: 
        :ref:`init`

    Reading data files within LAMMPS
        `LAMMPS read_data command  <https://docs.lammps.org/read_data.html>`_

    Specifying the units used in dimensions
        `LAMMPS units command  <https://docs.lammps.org/units.html>`_
