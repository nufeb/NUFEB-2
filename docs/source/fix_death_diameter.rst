.. index:: fix nufeb/death/diameter

fix nufeb/death/diameter command
================================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/death/diameter dead-group-ID diameter

* ID = the user-assigned name for the fix
* group-ID = the user-assigned group ID of atoms to which the fix is applied
* dead-group-ID = user-assigned group ID denoting dead atoms
* diameter = threshold diameter (meters) below which atoms are marked as dead

Examples
""""""""

.. code-block:: 

    group het type 1
    group dead type 2  
    
    fix death het nufeb/death/diameter dead 5e-7

Description
"""""""""""

Atoms below the threshold diameter are assigned to the *dead* group and no longer part of *group*.
