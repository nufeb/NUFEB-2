.. index:: fix nufeb/growth/aob

fix nufeb/growth/aob command
============================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/growth/aob nh4-ID nh4-Ks o2-ID o2-Ks no2-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* nh4-ID = substrate ID for ammonium
* nh4-Ks = half-velocity constant (Ks) for ammonium
* o2-ID = substrate ID for oxygen
* o2-Ks = half-velocity constant (Ks) for oxygen
* no2-ID = substrate ID for nitrite
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

   group aob type 1
   grid_style nufeb/chemostat 3 nh4 o2 no2 4e-6
   
   fix f_gaob aob nufeb/growth/aob nh4 1e-3 o2 5.4e-4 no2 growth 2.3727e-5 yield 0.15 maintain 1.505e-6 decay 1.27e-7
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as ammonia-oxidizing bacteria (AOB),
which have a spherical shape
(see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes that AOBs can can
perform ammonium oxidation to produce nitrite as part of their growth process.
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
  
  r1 & = \mu_{max} \cdot \frac{S_{nh4}}{S_{nh4} + Ks_{nh4}} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
  r2 & = b_{maint} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{nh4}, S_{o2}` are the local concentrations of ammonium and oxygen, respectively, at the grid cell in which atom resides
* :math:`Ks_{nh4}, Ks_{o2}` are the half-velocity constants for ammonium (*nh4-Ks*) and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient (*maintain*)
  
The new mass is then used to update the diameter of the atom.
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilisation (reaction) rates in all the affected grid cells:

.. math::
  
   R_{nh4} & = -\frac{1}{Y} \cdot r1 \cdot X
     
   R_{o2} & = -(\frac{3.42 - Y}{Y} \cdot r1 + r3) \cdot X
   
   R_{no2} & = \frac{3.42 - Y}{Y} \cdot r1 \cdot X
  
where:

* :math:`R_{nh4}, R_{o2}, R_{no2}` are the utilisation rates of ammonium, oxygen, and nitrite in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the AOB biomass density in grid cell

Restrictions
"""""""""""""
This fix is not compatible with the following command:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

----------

.. _ofiteru13:

**(Ofiteru, I.D., et al 2013)** Ofiteru, I.D., et al., Multi-scale modelling of bioreactor-separator system for wastewater
treatment with two-dimensional activated sludge floc dynamics, Water Research (2013)
   
  
  
  
  
  
  
  