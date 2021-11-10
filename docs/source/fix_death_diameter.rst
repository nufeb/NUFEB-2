.. index:: fix nufeb/death/diameter

fix nufeb/death/diameter command
================================

Syntax
""""""

.. parsed-literal::
    
     fix ID group nufeb/death/diameter dead-group diameter

* ID = the user-assigned name for the fix
* group = the user-assigned group of atoms/cells to which the fix is applied
* dead-group = user-assigned group denoting dead atoms/cells
* diameter = threshold diameter (meters) at which cells are marked as dead

Examples
""""""""

.. code-block:: 

    fix death het nufeb/death/diameter dead 5e-7

Description
"""""""""""

Cells below the threshold diameter are assigned to the *dead-dead* group and no longer part of *group*.
