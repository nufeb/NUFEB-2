.. index:: fix nufeb/growth/nob

fix nufeb/growth/nob command
==============================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/growth/nob o2-ID o2-Ks no2-ID no2-Ks no3-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* o2-ID = substrate ID for oxygen
* o2-Ks = half-velocity constant (Ks) for oxygen
* no2-ID = substrate ID for nitrite
* no2-Ks = half-velocity constant (Ks) for nitrite
* no3-ID = substrate ID for nitrate
* zero or more keyword/value pairs may be appended
* keyword = *growth* or *yield* or *decay* or *maintain* 

	.. parsed-literal::
	
	    *growth* value = maximum growth rate 
	    *yield* value = yield coefficient
	    *decay* value = decay rate
	    *maintain* value = maintenance coefficient

Examples
""""""""

.. code-block::

   #--- examples/biofilm-anammox ---#

   group nob type 1
   grid_style nufeb/chemostat 3 o2 no2 no3 4e-6
   
   fix f_gnob nob nufeb/growth/nob o2 6.8e-4 no2 1.3e-3 no3 growth 1.6782e-5 yield 0.041 maintain 0.694e-6 decay 1.27e-7
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as nitrite-oxidizing bacteria (NOB), which have a spherical shape
(see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes that NOBs can perform nitrite oxidation to produce nitrate
as part of their growth process.
Additionally, the model takes into account microbial decay and endogenous respiration processes.

The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The value of the substrate ID keyword *XX-ID* must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.
The following forward Euler method is implemented to update the mass 
(*m*) of each atom in the group:

.. math::

  m' & = m + \mu \cdot \Delta t
  
The specific growth rates :math:`\mu` is 
calculated based on the equations described in :ref:`(Ofiteru, I.D., et al, 2013) <ofiteru13>`: 

.. math::
  \mu & = r1 + r2 - b_{decay}
  
  r1 & = \mu_{max} \cdot \frac{S_{no2}}{S_{no2} + Ks_{no2}} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
  r2 & = b_{maint} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{no2}, S_{o2}` are the local concentrations of nitrite and oxygen, respectively, at the grid cell in which atom resides
* :math:`Ks_{no2}, Ks_{o2}` are the half-velocity constants for nitrite (*no2-Ks*) and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient (*maintain*)
  
The new mass is then used to update the diameter of the atom.
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilisation (reaction) rates in all the affected grid cells:

.. math::
  
   R_{no2} & = -\frac{1}{Y} \cdot r1 \cdot X
     
   R_{o2} & = -(\frac{1.15 - Y}{Y} r1 + r3) \cdot X
   
   R_{no3} & = \frac{1}{Y} \cdot r1 \cdot X
  
where:

* :math:`r_{no2}, r_{o2}, R_{no3}` are the utilisation rates of nitrite, oxygen, and nitrate in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the NOB biomass density in grid cell

Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

----------

.. _ofiteru13:

**(Ofiteru, I.D., et al 2013)** Ofiteru, I.D., et al., Multi-scale modelling of bioreactor-separator system for wastewater
treatment with two-dimensional activated sludge floc dynamics, Water Research (2013)
