.. index:: fix nufeb/growth/simple

fix nufeb/growth/simple command
===============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/simple sub-ID keyword value

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = ID of the substrate for atom growth
* keyword = *growth* 

	.. parsed-literal::
	
	    *growth* value = growth rate 

Examples
""""""""

.. code-block:: 

   group bac type 1
   grid_style nufeb/chemostat 1 glucose 0.001
   
   fix f_monod bac nufeb/growth/simple glucose growth 4e-4 

Description
"""""""""""

Perform exponential microbial growth (or decay) to the atoms defined in *group-ID*. 
The fix is called in each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The following forward Euler method is implemented to update the mass (*m*) of each atom in the group:

.. math::

  m'= m + \mu \Delta t
  
where :math:`\mu` is the growth rate of the atoms (*growth*). 
The new mass is then used to update other atom attributes. For a coccus-style atom,
its diameter changes accordingly. For a bacillus-style atom, the update is along
the length of the atom while the variations in its width (diameter) are negligible.