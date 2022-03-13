.. index:: fix nufeb/property/plasmid

fix nufeb/property/plasmid command
=====================================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/property/plasmid keyword value ...
    
* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* zero or more keyword/value pairs may be appended
* keyword = *max* or *init* or *dia* 

	.. parsed-literal::
	
	    *max* value = maximum # of plasmid
	    *init* value = # of initial plasmid
	    *dia* value = plasmid diameter

Examples
""""""""

.. code-block:: 

   fix p_plm ecoli nufeb/property/plasmid max 5 init 1 dia 1e-7
   
Description
"""""""""""

Assign plasmid attribute to each atom defined in *group-ID* .... TODO