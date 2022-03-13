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
Assign cell generation attribute to each atom defined in *group-ID* .... TODO