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
Assign cell cell cycle time attribute to each atom defined in *group-ID* .... TODO