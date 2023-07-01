.. index:: fix nufeb/growth/monod

fix nufeb/growth/monod command
===============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/monod sub-ID sub-Ks keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = ID of the substrate for atom growth
* sub-Ks = half-velocity constant (Ks) for the substrate
* zero or more keyword/value pairs may be appended
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
Impose a Monod-like growth process to the atoms defined in *group-ID*. The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The value of the substrate ID keyword *sub-ID* must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.
The following forward Euler method is implemented to update the mass (*m*) of each atom in the group:

.. math::

  m'= m + \mu \cdot \Delta t
  
The specific growth rate :math:`\mu` is calculated based on the Monod equation:

.. math::

  \mu = \mu_{max} \cdot \frac{S_{sub}}{S_{sub} + Ks_{sub}} - b_{decay}
  
where:

* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{sub}` is the local substrate concentration at the grid cell in which atom resides
* :math:`Ks_{sub}` is the half-velocity constant for the substrate (*sub-Ks*)
* :math:`b_{decay}` is the decay rate (*decay*)

The new mass is then used to update atom attributes. In the case of
:doc:`atom_style coccus <atom_vec_coccus>` is used,
the diameter changes accordingly.
For :doc:`atom_style bacillus <atom_vec_bacillus>`,
update affects the length of the bacilli.

If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilisation (reaction) rate R in all the affected grid cells:

.. math::

  R = \frac{1}{Y} \cdot \mu \cdot X
  
where:

* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the biomass density of the affected atoms in grid cell
