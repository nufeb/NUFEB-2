.. index:: nufeb/growth/anammox

fix nufeb/growth/anammox command
================================

Syntax
""""""

.. parsed-literal::

    fix ID group-ID nufeb/growth/anammox nh4-ID nh4-Ks o2-ID o2-Ks no2-ID no2-Ks no3-ID n2-ID keyword value ...
    
* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* nh4-ID = substrate ID for ammonium
* nh4-Ks = half-velocity constant (Ks) for ammonium
* o2-ID = substrate ID for oxygen
* o2-Ks = half-velocity constant (Ks) for oxygen
* no2-ID = substrate ID for nitrite
* no2-Ks = half-velocity constant (Ks) for nitrite
* no3-ID = substrate ID for nitrate
* n2-ID = substrate ID for nitrogen
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

   group ana type 1
   grid_style nufeb/chemostat 5 nh4 o2 no2 no3 n2 4e-6
   
   fix f_gana ana nufeb/growth/anammox nh4 7e-5 o2 1e-5 no2 5e-5 no3 n2 growth 9.26e-7 yield 0.159 maintain 3.5e-8 decay 3e-8
   
Description
""""""""""""""

Perform microbial growth to the atoms defined in *group-ID*.
The affected atoms are considered as anaerobic ammonia-oxidizing bacteria (anammox), which have a
spherical shape (see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes that anammox bacteria can transform ammonium and nitrite into nitrogen and low
nitrate production.
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
calculated based on the equations described in :ref:`(Hao, X., et al, 2002) <hao02>`: 

.. math::
  \mu & = r1 - r2 - b_{decay}
  
  r1 & = \mu_{max} \cdot \frac{S_{nh4}}{S_{nh4} + Ks_{nh4}} \cdot \frac{S_{no2}}{S_{no2} + Ks_{no2}} \cdot \frac{Ks_{o2}}{S_{o2} + Ks_{o2}}
  
  r2 & = b_{maint} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}

where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{nh4}, S_{no2}, S_{o2}` are the local concentrations of ammonium, nitrite, and oxygen, respectively, at the grid cell in which atom resides
* :math:`Ks_{nh4}, Ks_{no2}, Ks_{o2}` are the half-velocity constants for ammonium (*nh4-Ks*), nitrite (*no2-Ks*), and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient (*maintain*)
  

The new mass is then used to update the diameter.
If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also updates substrate utilisation (reaction) rates in all the affected grid cells:

.. math::
  
   R_{nh4} & = -\frac{1}{Y} \cdot r1 \cdot X
   
   R_{no2} & = -(\frac{1}{Y} + \frac{1}{1.14}) \cdot r1 \cdot X
   
   R_{no3} & = \frac{1}{1.14} \cdot r1 \cdot X

   R_{n2} & = \frac{2}{Y} \cdot r1 \cdot X
  
where:

* :math:`R_{nh4}, R_{no2}, R_{no3}, R_{n2}` are the utilisation rates of ammonium, nitrite, nitrate, and nitrogen in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the Anammox biomass density in grid cell


Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

----------

.. _hao02:

**(Hao, X., et al, 2002)** Hao, X., et al., 
Sensitivity analysis of a biofilm model describing a 
one-stage completely autotrophic nitrogen removal (CANON) process. Biotechnol. Bioeng (2002)