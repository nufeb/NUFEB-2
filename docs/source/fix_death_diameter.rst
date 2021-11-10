.. index:: fix nufeb/death/diameter

fix nufeb/death/diameter command
================================

Syntax
""""""

.. parsed-literal::
    
     fix ID group nufeb/death/diameter dead-group diameter

* ID = the user-assigned name for the fix
* group = the user-assigned group of atoms to which the fix is applied
* dead-group = user-assigned group denoting dead atoms
* diameter = threshold diameter (meters) below which atoms are marked as dead

Examples
""""""""

.. code-block:: 

    fix death het nufeb/death/diameter dead 5e-7

Description
"""""""""""

Atoms below the threshold diameter are assigned to the *dead-dead* group and no longer part of *group*.
