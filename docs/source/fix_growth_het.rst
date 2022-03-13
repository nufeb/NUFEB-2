.. index:: fix nufeb/growth/het

fix nufeb/growth/het command
============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/het sub-ID sub-Ks o2-ID o2-Ks no2-ID no2-Ks no3-ID no3-Ks keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = ID of the (organic) substrate for atom growth
* sub-Ks = half-velocity constant (Ks) of the substrate
* o2-ID = ID of the substrate (oxygen) for aerobic growth
* o2-Ks = half-velocity constant (Ks) of oxygen
* no2-ID = ID of the substrate (nitrite) for anaerobic growth
* no2-Ks = half-velocity constant (Ks) of nitrite
* no3-ID = ID of the substrate (nitrate) for anaerobic growth
* no3-Ks = half-velocity constant (Ks) of nitrate
* zero or more keyword/value pairs may be appended
* keyword = *growth* or *yield* or *decay* or *maintain* 

	.. parsed-literal::
	
	    *growth* value = maximum growth rate 
	    *yield* value = yield coefficient
	    *decay* value = decay rate
	    *maintain* value = maintenance coefficient
	    *epsyield* value = yield coefficient for EPS production 
	    *anoxic* value = reduction factor in anoxic condition
	    *epsdens* value = EPS density
         
Examples
""""""""

.. code-block:: 

   group het type 1
   grid_style nufeb/chemostat 4 sub o2 no2 no3 4e-6
   
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 growth 6.9e-5 yield 0.61
   
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 & 
   growth 6.9e-5 yield 0.61 epsyield 0.18 epsdens 30
      
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 & 
   growth 6.9e-5 yield 0.61 decay 9.2e-7 maintain 3.7e-6 epsyield 0.18 anoxic 0.6 epsdens 30
   
   
Description
"""""""""""
Perform microbial growth (or decay) to the atoms defined in *group-ID*. 
The affected atoms are considered as heterotrophic bacteria, 
in spherical shape
with outer mass and outer diameter for representing their EPS shells
(see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes heterotrophs grow by consuming organic substrate in 
oxygenated conditions or nitrate in anoxic denitrifying conditions,
and takes also into account the related decay and endogenous respiration processes.

The fix is called in each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The following forward Euler method is implemented to update the mass 
(*m*) and outer mass (*om*) of each atom in the group:

.. math::

  m' & = m + \mu \Delta t
  
  om' & = om + \mu_{EPS} \Delta t
  
The specific growth rates :math:`\mu` and EPS secretion rate :math:`\mu_{EPS}` are
calculated based on the equations described in :ref:`(Ofiteru, I.D., et al, 2013) <ofiteru13>`: 

.. math::
  \mu & = r1 + r2 + r3 - r4 - r5 - r6 - b_{decay}
  
  \mu_{EPS} & = \frac{Y_{EPS}}{Y} (r1 + r2 +r3)
    
  r1 & = \mu_{max} \frac{S_{sub}}{S_{sub} + Ks_{sub}} \frac{S_{o2}}{S_{o2} + Ks_{o2}} 
  
  r2 & = \eta \mu_{max} \frac{S_{sub}}{S_{sub} + Ks_{sub}} \frac{S_{no3}}{S_{no3} + Ks_{no3}} \frac{Ks_{o2}}{S_{o2} + Ks_{o2}} 
  
  r3 & = \eta \mu_{max} \frac{S_{sub}}{S_{sub} + Ks_{sub}} \frac{S_{no2}}{S_{no2} + Ks_{no2}} \frac{Ks_{o2}}{S_{o2} + Ks_{o2}} 
  
  r4 & = b_{maint} \frac{S_{o2}}{S_{o2} + Ks_{o2}} 
  
  r5 & = \frac{1}{1.17} \eta b_{maint} \frac{S_{no2}}{S_{no2} + Ks_{no2}} \frac{Ks_{o2}}{S_{o2} + Ks_{o2}} 
  
  r6 & = \frac{1}{2.86} \eta b_{maint} \frac{S_{no3}}{S_{no3} + Ks_{no3}} \frac{Ks_{o2}}{S_{o2} + Ks_{o2}} 
  
where:

* :math:`b_{decay}` is the decay rate of the atoms (*decay*)
* :math:`Y` is the yield coefficient of the atoms (*yield*)
* :math:`Y_{EPS}` is the yield coefficient for EPS secretion of the atoms (*epsyield*)
* :math:`\mu_{max}` is the maximum growth rate of the atoms (*growth*)
* :math:`S_{sub}, S_{o2}, S_{no2}, S_{no3}` are the local concentrations of organic substrate, oxygen, nitrite and nitrate at the grid cell in which atom resides, respectively
* :math:`Ks_{sub}, Ks_{o2}, Ks_{no2}, Ks_{no3}` are the half-velocity constants of the substrate (*sub-Ks*), oxygen (*o2-Ks*), nitrite (*no2-Ks*) and nitrate (*no3-Ks*), respectively
* :math:`\eta` is the reduction factor of the atoms in anoxic condition (*anoxic*)
* :math:`b_{maint}` is the maintenance coefficient of the atoms (*maintain*)

The new mass and outer mass are then used to update diameter and outer diameter of the atom, respectively. 
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilization (reaction) rates in all the affected grid cells. 
The rates are related to the specific growth rate and yield as follows:

.. math::

  r_{sub} & = -\frac{1}{Y} (r1 + r2 + r3) X
  
  r_{o2} & = -(\frac{1-Y-Y_{EPS}}{Y} r1 + r4) X
  
  r_{no3} & = -(\frac{1-Y-Y_{EPS}}{2.86 Y} r2  + r5) X
    
  r_{no2} & = -(\frac{1-Y-Y_{EPS}}{1.17 Y} r3  + r6) X
  
  
where:

* :math:`r_{sub}, r_{o2}, r_{no2}, r_{no3}` are the utilization rates of organic substrate, oxygen, nitrite and nitrate in the affected grid cells, respectively
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

