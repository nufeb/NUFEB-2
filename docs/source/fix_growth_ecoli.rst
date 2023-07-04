.. index:: fix nufeb/growth/ecoli

fix nufeb/growth/ecoli command
====================================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/growth/ecoli suc-ID suc-Ks o2-ID o2-Ks co2-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* suc-ID = substrate ID for sucrose
* suc-Ks = half-velocity constant (Ks) for sucrose
* o2-ID = substrate ID for oxygen
* o2-Ks = half-velocity constant (Ks) for oxygen
* co2-ID = substrate ID for carbon dioxide
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

   #--- examples/Pub-J.Sakkos-2023-Phototroph ---#

   group ecoli type 1
   grid_style nufeb/chemostat 3 suc o2 co2 4e-6
   
   fix f_gecoli aob nufeb/growth/ecoli suc 3.6 o2 0.001 co2 growth 2.7e-4 yield 0.43 maintain 9.5e-7 decay 2e-5
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as *Escherichia coli* (*E.coli*).
The model assumes that *E.coli* can utilise sucrose and oxygen to produce co2.
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
calculated based on the equations described in :ref:`(Sakkos, J., et al, 2023) <sakkos23>`:

.. math::
  \mu & = r1 - r2 - b_{decay}

  r1 & = \mu_{max} \cdot \frac{S_{suc}}{S_{suc} + Ks_{suc}} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}

  r2 & = b_{maint} \cdot \frac{S_{o2}}{S_{o2} + Ks_{o2}}

where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{suc}, S_{o2}` are the local concentrations of sucrose and oxygen, respectively, at the grid cell in which atom resides
* :math:`Ks_{suc}, Ks_{o2}` are the half-velocity constants for sucrose (*suc-Ks*) and oxygen (*o2-Ks*), respectively
* :math:`b_{maint}` is the maintenance coefficient (*maintain*)

The new mass is then used to update atom attributes. In the case of
:doc:`atom_style coccus <atom_vec_coccus>` is used,
the diameter changes accordingly.
For :doc:`atom_style bacillus <atom_vec_bacillus>`,
update affects the length of the bacilli.

If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also updates substrate utilisation (reaction) rates in all the affected grid cells:

.. math::

   R_{suc} & = -\frac{1}{Y} \cdot r1 \cdot X

   R_{o2} & = -0.399  \cdot (r1 + r2) \cdot X

   R_{co2} & = -0.2  \cdot (r1 + r2) \cdot X


where:

* :math:`R_{suc}, R_{o2}, R_{co2}` are the utilisation rates of sucrose, oxygen, carbon dioxide in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the *E.coli* biomass density in grid cell

----------

.. _sakkos23:

**(Sakkos, J., et al, 2023)** Sakkos, J., et al.,
Predicting partner fitness based on spatial structuring in a light-driven microbial community.
PLoS Comput. Biol. (2023)