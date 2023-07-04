.. index:: fix nufeb/growth/cyano

fix nufeb/growth/cyano command
===============================

Syntax
""""""

.. parsed-literal::

    fix ID group-ID nufeb/growth/cyano light-ID light-Ks o2-ID co2-ID co2-Ks suc-ID gco2-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* suc-ID = substrate ID for light
* suc-Ks = half-velocity constant (Ks) for light
* o2-ID = substrate ID for oxygen
* co2-ID = substrate ID for carbon dioxide
* co2-Ks = half-velocity constant (Ks) for carbon dioxide
* suc-ID = substrate ID for sucrose
* gco2-ID = substrate ID for gaseous carbon dioxide
* zero or more keyword/value pairs may be appended
* keyword = *growth* or *yield* or *decay* or *maintain*

	.. parsed-literal::

	    *growth* value = maximum growth rate
	    *yield* value = yield coefficient
	    *decay* value = decay rate
	    *maintain* value = maintenance coefficient
	    *suc_exp* value = sucrose export rate


Examples
""""""""

.. code-block::

   #--- examples/Pub-J.Sakkos-2023-Phototroph ---#

   group ecoli type 1
   grid_style nufeb/chemostat 3 suc o2 co2 4e-6

   fix f_gcyano aob nufeb/growth/cyano light 3.5e-04 o2 co2 2e-4 suc gco2 &
   growth 1.67e-05 yield 0.55 suc_exp 0.285

Description
""""""""""""""

Perform microbial growth to the atoms defined in *group-ID*.
The affects atoms are considered as *Synechococcus elongatus PCC 7942* -
an engineered cyanobacterial strain that can secrete sucrose by utilising light and carbon.

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
  \mu & = r1 \cdot r2 - b_{decay} - b_{maint}

  r1 & = \mu_{max} \cdot \frac{S_{light}}{S_{light} + Ks_{light}} \cdot \frac{S_{co2}}{S_{co2} + Ks_{co2}}

  r2 & = 0.141 \cdot e^{\frac{-suc\_exp}{0.063}} + 0.9

where:

* :math:`b_{decay}` is the decay rate (*decay*)
* :math:`b_{maint}` is the maintenance rate (*maintain*)
* :math:`\mu_{max}` is the maximum growth rate (*growth*)
* :math:`S_{light}, S_{co2}` are the local concentrations of light and carbon dioxide, respectively, at the grid cell in which atom resides
* :math:`Ks_{light}, Ks_{co2}` are the half-velocity constants for light (*light-Ks*) and carbon dioxide (*co2-Ks*), respectively
* :math:`r2`  is an empirical fit for growth reduction with respect to IPTG induction of the sucrose secretion machinery

The new mass is then used to update atom attributes. In the case of
:doc:`atom_style coccus <atom_vec_coccus>` is used,
the diameter changes accordingly.
For :doc:`atom_style bacillus <atom_vec_bacillus>`,
update affects the length of the bacilli.

If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also updates substrate utilisation (reaction) rates in all the affected grid cells:

.. math::

   \psi & = r1 \cdot (-3.4897 \cdot e^{\frac{-suc\_exp}{0.048}} + 3.4092)

   R_{light} & = -\frac{1}{Y} \cdot (r1 + \psi) \cdot X

   R_{co2} & = -\frac{1}{Y} \cdot (r1 + \psi) \cdot X

   R_{o2} & = \frac{0.727}{Y} \cdot (r1 + \psi) \cdot X - 0.1 \cdot b_{maint} \cdot X

   R_{suc} & = \frac{0.65}{Y} \cdot (r1 + \psi) \cdot X

where:

* :math:`\psi` is the metabolic flux due to sucrose secretion
* :math:`R_{light}, R_{co2}, R_{o2}, R_{suc}` are the utilisation rates of sucrose, carbon dioxide, oxygen, and sucrose in the affected grid cells, respectively
* :math:`Y` is the yield coefficient (*yield*)
* :math:`X` is the *E.coli* biomass density in grid cell

----------

.. _sakkos23:

**(Sakkos, J., et al, 2023)** Sakkos, J., et al.,
Predicting partner fitness based on spatial structuring in a light-driven microbial community.
PLoS Comput. Biol. (2023)