.. index:: fix nufeb/growth/aob

fix nufeb/growth/aob command
============================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/growth/aob nh4-ID nh4-Ks o2-ID o2-Ks no2-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* nh4-ID = ID of the substrate (ammonium) for atom growth 
* nh4-Ks = half-velocity constant (Ks) for ammonium
* o2-ID = ID of the substrate (oxygen)
* o2-Ks = half-velocity constant (Ks) for oxygen
* no2-ID = ID of the substrate (nitrite)
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

   group aob type 1
   grid_style nufeb/chemostat 3 nh4 o2 no2 4e-6
   
   fix f_gaob aob nufeb/growth/aob nh4 1e-3 o2 5.4e-4 no2 growth 2.3727e-5 yield 0.15 maintain 1.505e-6 decay 1.27e-7
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as ammonia-oxidizing bacteria (AOB), 
in spherical shape without outer mass and outer diameter
(see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes AOBs can oxidize ammonium to nitrite for their growth,
and takes also into account the related decay and endogenous respiration processes.

The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The following forward Euler method is implemented to update the mass 
(*m*) of each atom in the group:

.. math::

  m' & = m + \mu \Delta t
  
The specific growth rates :math:`\mu` is 
calculated based on the equations described in :ref:`(Ofiteru, I.D., et al, 2013) <ofiteru13>`: 

.. math::
  \mu & = r1 + r2 - b_{decay}
  
  r1 & = \mu_{max} \frac{S_{nh4}}{S_{nh4} + Ks_{nh4}} \frac{S_{o2}}{S_{o2} + Ks_{o2}} 
  
  r2 & = b_{maint} \frac{S_{o2}}{S_{o2} + Ks_{o2}} 
  
where:

* :math:`b_{decay}` is the decay rate of the atoms (*decay*)
* :math:`\mu_{max}` is the maximum growth rate of the atoms (*growth*)
* :math:`S_{nh4}, S_{o2}` are the local concentrations of ammonium and oxygen at the grid cell in which atom resides, respectively
* :math:`Ks_{nh4}, Ks_{o2}` are the half-velocity constants for ammonium (*nh4-Ks*) and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient of the atoms (*maintain*)
  
The new mass is then used to update diameter of the atom. 
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilization (reaction) rates in all the affected grid cells. 
The rates are related to the specific growth rate and yield as follows:

.. math::
  
   r_{nh4} & = -\frac{1}{Y} r1 X
     
   r_{o2} & = -(\frac{3.42 - Y}{Y} r1 + r3) X
   
   r_{no2} & = \frac{3.42 - Y}{Y} r1 X
  
where:

* :math:`r_{nh4}, r_{o2}, r_{no2}` are the utilization rates of ammonium, oxygen, and nitrite in the affected grid cells, respectively
* :math:`X` is the biomass density in grid cell 

Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

* :doc:`grid_style simple <grid_style_simple>`

----------

.. _ofiteru13:

**(Ofiteru, I.D., et al 2013)** Ofiteru, I.D., et al., Multi-scale modelling of bioreactor-separator system for wastewater
treatment with two-dimensional activated sludge floc dynamics, Water Research (2013)
   
  
  
  
  
  
  
  