.. index:: fix nufeb/growth/het

fix nufeb/growth/het command
============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/het sub-ID sub-Ks o2-ID o2-Ks no2-ID no2-Ks no3-ID no3-Ks keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = substrate ID for organic substrate
* sub-Ks = half-velocity constant (Ks) for organic substrate
* o2-ID = substrate ID for oxygen
* o2-Ks = half-velocity constant (Ks) for oxygen
* no2-ID = substrate ID for nitrite
* no2-Ks = half-velocity constant (Ks) for nitrite
* no3-ID = substrate ID for nitrate
* no3-Ks = half-velocity constant (Ks) for nitrate
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

   #--- examples/biofilm-heterotroph ---#

   group het type 1
   grid_style nufeb/chemostat 4 sub o2 no2 no3 4e-6
   
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 growth 6.9e-5 yield 0.61
   
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 & 
   growth 6.9e-5 yield 0.61 epsyield 0.18 epsdens 30
      
   fix f_ghet het nufeb/growth/het sub 4e-3 o2 2e-4 no2 3e-4 no3 3e-4 & 
   growth 6.9e-5 yield 0.61 decay 9.2e-7 maintain 3.7e-6 epsyield 0.18 anoxic 0.6 epsdens 30
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*.
The affected atoms are considered as heterotrophic bacteria, having a spherical shape (:doc:`atom_style coccus <atom_vec_coccus>`)
with outer mass and outer diameter to represent their EPS (Extracellular Polymeric Substance) shells.
The model assumes that heterotrophs grow by consuming organic substrate in
oxygenated conditions or nitrate in anoxic denitrifying conditions.
Additionally, the model takes into account microbial decay and endogenous respiration processes.

The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.
The value of the substrate ID keyword *XX-ID* must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.
The following forward Euler method is implemented to update the mass 
(*m*) and outer mass (*om*) of each atom in the group:

.. math::

  m' & = m + \mu \cdot \Delta t
  
  om' & = om + \mu_{EPS} \cdot \Delta t
  
The specific growth rates :math:`\mu` and EPS secretion rate :math:`\mu_{EPS}` are
calculated based on the equations described in :ref:`(Ofiteru, I.D., et al, 2013) <ofiteru13>`:

.. math::

  \mu & = r1 + r2 + r3 - r4 - r5 - r6 - b_{decay}

  \mu_{EPS} & = \frac{Y_{EPS}}{Y} \cdot (r1 + r2 +r3)

  r1 & = \mu_{max} \cdot \frac{S_{sub}}{S_{sub} + Ks_{sub}} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
  r2 & = \eta \mu_{max} \cdot \frac{S_{sub}}{S_{sub} + Ks_{sub}} \cdot \frac{S_{no3}}{S_{no3} + Ks_{no3}} \cdot \frac{Ks_{o2}}{S_{o2} + Ks_{o2}}
  
  r3 & = \eta \mu_{max} \cdot \frac{S_{sub}}{S_{sub} + Ks_{sub}} \cdot \frac{S_{no2}}{S_{no2} + Ks_{no2}} \cdot \frac{Ks_{o2}}{S_{o2} + Ks_{o2}}
  
  r4 & = b_{maint} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}
  
  r5 & = \frac{1}{1.17} \cdot \eta \cdot b_{maint} \cdot \frac{S_{no2}}{S_{no2} + Ks_{no2}} \cdot \frac{Ks_{o2}}{S_{o2} + Ks_{o2}}
  
  r6 & = \frac{1}{2.86} \cdot \eta \cdot b_{maint} \cdot \frac{S_{no3}}{S_{no3} + Ks_{no3}} \cdot \frac{Ks_{o2}}{S_{o2} + Ks_{o2}}
  
where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`Y` is the yield coefficient (*yield*)
* :math:`Y_{EPS}` is the yield coefficient for EPS secretion (*epsyield*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{sub}, S_{o2}, S_{no2}, S_{no3}` are the local concentrations of organic substrate, oxygen, nitrite and nitrate, respectively, at the grid cell in which atom resides
* :math:`Ks_{sub}, Ks_{o2}, Ks_{no2}, Ks_{no3}` are the half-velocity constants for organic substrate (*sub-Ks*), oxygen (*o2-Ks*), nitrite (*no2-Ks*) and nitrate (*no3-Ks*), respectively
* :math:`\eta` is the reduction factor of the atoms in anoxic condition (*anoxic*)
* :math:`b_{maint}` is the maintenance coefficient (*maintain*)

The new mass and outer mass are then used to update the diameter and outer diameter of the atoms.
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilisation (reaction) rates in all the affected grid cells:

.. math::

  R_{sub} & = -\frac{1}{Y} \cdot (r1 + r2 + r3) \cdot X
  
  R_{o2} & = -(\frac{1-Y-Y_{EPS}}{Y} \cdot r1 + r4) \cdot X
  
  R_{no3} & = -(\frac{1-Y-Y_{EPS}}{2.86 Y} \cdot r2  + r5) \cdot X
    
  R_{no2} & = -(\frac{1-Y-Y_{EPS}}{1.17 Y} \cdot r3  + r6) \cdot X
  
  
where:

* :math:`R_{sub}, R_{o2}, R_{no2}, R_{no3}` are the utilisation rates of organic substrate, oxygen, nitrite and nitrate in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the heterotrophs biomass density in grid cell


Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

----------

.. _ofiteru13:

**(Ofiteru, I.D., et al 2013)** Ofiteru, I.D., et al., Multi-scale modelling of bioreactor-separator system for wastewater
treatment with two-dimensional activated sludge floc dynamics, Water Research (2013)

