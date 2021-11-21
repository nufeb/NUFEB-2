.. index:: fix nufeb/growth/monod

fix nufeb/growth/monod command
===============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/monod sub-ID Ks keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = ID of the substrate which is utilized by the group of atoms
* sub-Ks = half-velocity constant (Ks) of the substrate
* keyword = *growth* or *yield* or *decay* 

	.. parsed-literal::
	
	    *growth* value = maximum growth rate 
	    *yield* value = yield coefficient
	    *decay* value = decay rate

         
Examples
""""""""

.. code-block:: 

   group bac type 1
   grid_style nufeb/chemostat 1 glucose 0.001
   
   fix f_monod bac nufeb/growth/monod glucose 3.5e-5 growth 4e-4 yield 0.61 decay 2e-5 maintain 1e-5
   fix f_monod bac nufeb/growth/monod glucose 3.5e-5 growth 4e-4 yield 0.61 

Description
"""""""""""
Perfrom microbial growth (or decay) to the atoms defined in *group-ID*. The fix is called in each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid properties.
The following forward Euler method is implemented to update the mass of each atom in the group:

.. math::

  m'= m + \mu \Delta t
  
The specific growth rate :math:`\mu` is calcualted based on the Monod equation:

.. math::

  \mu = \mu_{max} \frac{S_{sub}}{S_{sub} + Ks_{sub}} - b_{decay}
  
where:

* :math:`\mu_{max}` is the maximum growth rate of the atoms (*growth*)
* :math:`S_{sub}` is the local substrate concentration at the grid cell in which atom resides
* :math:`Ks_{sub}` is the half-velocity constant of the substrate (*sub-Ks*)
* :math:`b_{decay}` is the decay rate of the atoms (*decay*)

The new mass is then used to update other atom properties. For a coccus-style atom,
its diameter changes accordingly. For a bacillus-style atom, the update is along
the length of the atom while the variations in its width (diameter) are neglible.

If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilization (reaction) rate for all the affected grid cells. 
The rate is related to the specific growth rate as follows:

.. math::

  r = \frac{1}{Y} \mu X
  
where:

* :math:`Y` is the yield coefficient of the atoms (*yield*)
* :math:`X` is the biomass density in grid cell 
