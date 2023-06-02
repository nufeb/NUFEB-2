.. index:: fix nufeb/growth/energy

fix nufeb/growth/energy command
===============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/growth/energy file keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* file = name of data file to read in
* zero or more keyword/value pairs may be appended
* keyword = *temperature* or *alfa* or *beta* or *mw_biomass*

	.. parsed-literal::
	
	    *temperature* value = temperature (default: 298.15K)
	    *alfa* value = relaxation parameter for maintenance (default: 1.2)
	    *beta* value = relaxation parameter for local environment (default: 0.8)
	    *mw_biomass* value = molecular weight of biomass (default: 24.6g/mol)

         
Examples
""""""""

.. code-block::

   group bac type 1

   fix f_bac bac nufeb/growth/energy energy.in
   fix f_bac bac nufeb/growth/energy energy.in alfa 1.1 beta 0.9

Description
"""""""""""
Perform energy-based microbial growth to the atoms defined in *group-ID*.
The fix implements the model presented in
:ref:`(Gogulancea, V., et al 2019) <Gogulancea19>`,
and is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.

The growth process is modeled using thermodynamic
principles, enabling the estimation of growth yields according
to the chemical energy of the environment. More specifically,
the growth yield *Y* (mol-X/mol-eD) is estimated as:

.. math::
  Y & = \frac{\Delta G_{cat}{\Delta G_{ana} + \Delta G_{dis}}

where:

* :math:`\Delta G_{cat}` is the free energy required for microbial anabolic pathway
* :math:`\Delta G_{ana}` is the energy available from microbial catabolic pathway
* :math:`\Delta G_{dis}` is the *Dissipation Energy* for maintenance requirements defined in the data *file*

The anabolic and catabolic energies are calculated using the reaction quotient
and the standard *Substrate Gibbs Energy* change (defined in the data *file*) at certain
*temperature*.

In addition, the metabolic rate *q_cat* and maintenance energy *m_req* of each
microbial group are calculated as:

.. math::

  q_cat & = q_max * \frac{S_{sub-X}}{S_{sub-X} + Ks_{sub-X}}

  m_req & = \frac{0.00125}{\Delta G_{cat}}

*q_cat* depends on the maximum substrate *Uptake Rate* defined in the data *file*
and the availability of the limiting substrate. A constant value of 0.00125KJ/mol.s is assumed
for the maintenance energy of all the microbial groups.



Given the above, the growth model assumes mixed kinetic–thermodynamic
limitation, considering three possible scenarios for
microbial growth:




----------

Format of data file
"""""""""""""""""""""

A data file has a header and a body. The header appears first. The first line of the
header is always skipped; it typically contains a description of the file.
Lines can have a trailing comment starting with ‘#’ that is ignored. If the line is blank (only whitespace
after comment is deleted), it is skipped.

The following header is required:

* *groups* = # of groups to apply the fix to

example:

  .. parsed-literal::

    Additional parameters used in fix nufeb/growth/energy

     2 groups

The body of the file contains 9 sections.
The first line of a section has only
a keyword. The next line is skipped. The remaining lines of the section contain values.
The number of lines depends on the section keyword as described below. Zero or more
blank lines can be used between sections. Sections can appear in any order.

These are the section keywords for the body of the file.

* *Uptake Rate, Decay Rate, Substrate Gibbs Energy, Biomass Gibbs Energy, Dissipation Energy*
* *Ks Coeffs, Catabolic Coeffs, Anabolic Coeffs, Decay Coeffs*


----------

*Uptake Rate* section:

* one line per group
* line syntax: group-ID value
* unit: mol-eD/mol-X·s （eD = electron donor, X = biomass）

* example:

  .. parsed-literal::

       Uptake Rate

            AOB 5.7175e-5
            NOB 1.1078e-4

Define maximum substrate uptake rate of each group.
The value in this section must be greater than or equal to 0.

----------

*Decay Rate* section:

* one line per group
* line syntax: group-ID value
* unit: s\^-1

* example:

  .. parsed-literal::

       Decay Rate

            AOB 2.778e-6
            NOB 2.444e-6

Define decay rate of each group.
The value in this section must be greater than or equal to 0.

----------

*Substrate Gibbs Energy* section:

* one line per substrate
* line syntax: substrate-ID value
* unit: KJ/mol

* example:

  .. parsed-literal::

     Substrate Gibbs Energy

         nh3 -26.57
         no2 -32.20
         no3 -103.70
         o2  16.4
         co2 -586.7
         h2o -237.18

Assign substrate Gibbs free energy to each substrate defined in the
:doc:`grid_style chemostat <grid_style_chemostat>`
or :doc:`grid_style simple <grid_style_simple>` command.
Substrate lines can come in any order.

----------

*Biomass Gibbs Energy* section:

* one line per group
* line syntax: group-ID value
* unit: KJ/mol

* example:

  .. parsed-literal::

     Biomass Gibbs Energy

         AOB -67
         NOB -67

Define biomass Gibbs free energy of each group.

----------

*Dissipation Energy* section:

* one line per group
* line syntax: group-ID value
* unit: KJ/mol

* example:

  .. parsed-literal::

     Dissipation Energy

         AOB -3500
         NOB -3500

Define dissipation energy of each group. The value indicates the amount of energy
dissipated for microbial maintenance requirements.

----------

*Ks Coeffs* section:

* one line per group
* line syntax: group-ID sub-1 value sub-2 value ... sub-N value
* unit: kg/m\^3

* example:

  .. parsed-literal::

     Ks Coeffs

         AOB nh3 3.6e-5   o2 3e-5
         NOB no2 1.81e-7  o2 6.02e-5

Define half-velocity coefficients (Ks) of each group.
*sub-X* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`
or :doc:`grid_style simple <grid_style_simple>` command.
The value in this section must be positive.

----------

*Catabolic Coeffs* section:

* one line per group
* line syntax: group-ID sub-1 value sub-2 value ... sub-N value
* unit: kg/m\^3

* example:

  .. parsed-literal::

     Catabolic Coeffs

         AOB nh3 -1  no2 1  o2 -1.5  h2o 1  h 1
         NOB no2 -1  no3 1  o2 -0.5

Define microbial catabolic coefficients of each group.
The coefficients indicate the stoichiometric relationship between the
substrates and products in the microbial catabolic reaction.
*sub-X* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`
or :doc:`grid_style simple <grid_style_simple>` command.


----------

*Anabolic Coeffs* section:

* one line per group
* line syntax: group-ID sub-1 value sub-2 value ... sub-N value
* unit: kg/m\^3

* example:

  .. parsed-literal::

     Anabolic Coeffs

         AOB nh3 -0.9  no2 0.7  co2 -1  h2o 1.1  h -1
         NOB no2 -2.9  no3 2.7  co2 -1  h2o 0.2  h -1

Define microbial anabolic coefficients of each group.
The coefficients indicate the stoichiometric relationship between the
substrates and products in the microbial anabolic reaction.
*sub-X* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`
or :doc:`grid_style simple <grid_style_simple>` command.

----------

*Decay Coeffs* section:

* one line per group
* line syntax: group-ID sub-1 value sub-2 value ... sub-N value
* unit: kg/m\^3

* example:

  .. parsed-literal::

     Decay Coeffs

         AOB nh3 0.2   co2 1
         NOB nh3 0.2   co2 1

Define microbial decay coefficients of each group.
The coefficients indicate the relative amount of substrates released to the environment
during the microbial decay.
*sub-X* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`
or :doc:`grid_style simple <grid_style_simple>` command.

----------

.. _Gogulancea19:

**(Gogulancea, V., et al 2019)** Gogulancea, V., et al.,
Individual Based Model Links Thermodynamics, Chemical Speciation and
Environmental Conditions to Microbial Growth, Frontiers in Microbiology (2019)