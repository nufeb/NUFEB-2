.. index:: nufeb/growth/anammox

fix nufeb/growth/anammox command
================================

Syntax
""""""

.. parsed-literal::

    fix ID group-ID nufeb/growth/anammox nh4-ID nh4-Ks o2-ID o2-Ks no2-ID no2-Ks no3-ID keyword value ...
    
* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* nh4-ID = ID of the substrate (ammonium) for atom growth 
* nh4-Ks = half-velocity constant (Ks) of ammonium
* o2-ID = ID of the substrate (oxygen)
* o2-Ks = half-velocity constant (Ks) of oxygen
* no2-ID = ID of the substrate (nitrite)
* no2-Ks = half-velocity constant (Ks) of nitrite
* no3-ID = ID of the substrate (nitrate)
* keyword = *growth* or *yield* or *decay* or *maintain* 

	.. parsed-literal::
	
	    *growth* value = maximum growth rate 
	    *yield* value = yield coefficient
	    *decay* value = decay rate
	    *maintain* value = maintenance coefficient

Examples
""""""""

.. code-block:: 

   group ana type 1
   grid_style nufeb/chemostat 4 nh4 o2 no2 no3 4e-6
   
   fix f_gana ana nufeb/growth/anammox nh4 7e-5 o2 1e-5 no2 5e-5 no3 growth 9.26e-7 yield 0.159 maintain 3.5e-8 decay 3e-8
   
 Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as anaerobic ammonia-oxidizing bacteria (anammox), 
in spherical shape without outer mass and outer diameter
(see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes anammox can transform ammonium and nitrite into nitrigen gas as well as low
nitrate production, and it takes also into account the related decay and endogenous respiration processes.

The fix is called in each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid properties.
The following forward Euler method is implemented to update the mass 
(*m*) of each atom in the group:
 
  m' & = m + \mu \Delta t
  
The specific growth rates :math:`\mu` is 
calculated based on the equations described in :ref:`(Hao, X., et al, 2002) <hao02>`: 

.. math::
  \mu & = r1 + r2 - b_{decay}
  
  r1 & = \mu_{max} \frac{S_{nh4}}{S_{nh4} + Ks_{nh4}} \frac{S_{no2}}{S_{no2} + Ks_{no2}} \frac{Ks_{o2}}{S_{o2} + Ks_{o2}} 
  
  r2 & = b_{maint} \frac{S_{o2}}{S_{o2} + Ks_{o2}} 

where:

* :math:`b_{decay}` is the decay rate of the atoms (*decay*)
* :math:`\mu_{max}` is the maximum growth rate of the atoms (*growth*)
* :math:`S_{nh4}, S_{no2}, S_{o2}` are the local concentrations of ammonium, nitrite, and oxygen at the grid cell in which atom resides, respectively
* :math:`Ks_{nh4}, Ks_{no2}, Ks_{o2}` are the half-velocity constants of ammonium (*nh4-Ks*), nitrite (*no2-Ks*), and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient of the atoms (*maintain*)
  

The new mass is then used to update diameter of the atom. 
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilization (reaction) rates in all the affected grid cells. 
The rates are related to the specific growth rate and yield as follows:

.. math::
  
   r_{nh4} & = -\frac{1}{Y} r1 X
   
   r_{no2} & = -(\frac{1}{Y} + \frac{1}{1.14}) r1 X
   
   r_{no3} & = \frac{1}{1.14} r1 X
  
where:

* :math:`r_{nh4}, r_{no2}, r_{no3}` are the utilization rates of ammonium, nitrite, and nitrate in the affected grid cells, respectively
* :math:`X` is the biomass density in grid cell 


Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

* :doc:`grid_style simple <grid_style_simple>`

----------

.. _hao02:

**(Hao, X., et al, 2002)** Hao, X., et al., 
Sensitivity analysis of a biofilm model describing a 
one-stage completely autotrophic nitrogen removal (CANON) process. Biotechnol. Bioeng (2002)