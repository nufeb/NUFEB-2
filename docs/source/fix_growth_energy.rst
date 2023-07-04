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
	    *alfa* value = maximum relaxation parameter for maintenance (default: 1.2)
	    *beta* value = minimum relaxation parameter for maintenance (default: 0.8)
	    *mw_biomass* value = molecular weight of biomass (default: 24.6g/mol)

         
Examples
""""""""

.. code-block::

   #--- examples/Pub-V.Gogulancea-2019-IbM-Thermo ---#

   group AOB type 1
   group NOB type 2

   fix f_gaob AOB nufeb/growth/energy energy.in
   fix f_gnob NOB nufeb/growth/energy energy.in

Description
"""""""""""
Perform energy-based microbial growth to the atoms defined in *group-ID*.
The fix implements the model presented in
:ref:`(Gogulancea, V., et al 2019) <Gogulancea19>` :ref:`(González-Cabaleiro, R., et al 2015) <Cabaleiro15>`.
It is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.

In contract to the other fix nufeb/growth/\*, where the growth model of
each microbial functional group is directly embedded into the source code,
the fix employs a generalisable framework allows users to create model with
customisable parameters without modifying the underlying code.

The fix requires a separate data *file* containing information about kinetic and thermodynamic
parameters. The file must be in ASCII text format,
and the specific format of the data file is detailed in the next section.

.. note::

   This command is capable of automatic unit conversion as
   some parameters used are in moles (mol).
   However, the molecular weight of each substrate must be provided in :doc:`grid_modify <grid_modify>`
   using the *mw* keyword.

The growth process is modeled using thermodynamic
principles, enabling the estimation of growth yields  :math:`Y` according
to the chemical energy of the environment:

.. math::

  Y & = \frac{ \Delta G_{cat}}{ \Delta G_{ana} + \Delta G_{dis} }

where:

* :math:`\Delta G_{cat}` is the free energy required for microbial anabolic pathway, using its absolute value
* :math:`\Delta G_{ana}` is the energy available from microbial catabolic pathway
* :math:`\Delta G_{dis}` is the *Dissipation Energy* for maintenance requirements defined in the data *file*

:math:`\Delta G_{cat}` and :math:`\Delta G_{ana}` are calculated using
gas constant :math:`R` (0.0083144KJ/mol.K), *temperature* :math:`T`, reaction quotient :math:`Q`,
as well as the standard *Substrate Gibbs Energy* change defined in the data *file*:

.. math::

    \Delta G & = \Delta G^{o} + R \cdot T \cdot ln(Q)

In addition, the catabolic rate :math:`q_{cat}` (mol-eD/mol-X·s),
maintenance rate :math:`m_{req}` (mol-eD / mol-X·s),
and specific growth rate :math:`\mu` of the
functional group are calculated as:

.. math::

  & q_{cat} =  q_{max} \cdot \frac{S_{sub_i}}{S_{sub_i} + Ks_{sub_i}}

  & m_{req} = -\frac{m_{G}}{\Delta G_{cat}}

  & \mu  = (q_{cat} - m_{req}) \cdot Y

where:

* :math:`q_{cat}` depends on the maximum substrate *Uptake Rate* (:math:`q_{max}`) defined in the data *file* and the availability of the limiting substrate concentrations. Monod-like microbial kinetic model is adopted.
* :math:`m_{req}` is a fraction of the catabolic energy which is diverted from growth to maintenance purposes. It is assumed only temperature dependent and a constant value :math:`m_{G}` = 0.00125KJ/mol.s is considered.
* :math:`\mu` is the positive net growth occurring when the rate of energy harvest exceeds that required for the maintenance.

Given the above results, the growth model assumes mixed kinetic–thermodynamic
limitation, considering three possible scenarios for updating the biomass :math:`m` of each atom in the group.

.. math::

    \frac{dm}{dt} & = \mu  \cdot m   &  \text{if } q_{cat} > \alpha \cdot m_{req}

    \frac{dm}{dt} & =  0   &  \text{if } \beta \cdot m_{req} \le q_{cat} \le \alpha \cdot m_{req}

    \frac{dm}{dt} & = -D_{decay} \cdot \frac{(m_{req} - q_{cat})}{m_{req}} \cdot Y \cdot m   &  \text{if } q_{cat} < \beta \cdot m_{req}

If :doc:`fix nufeb/diffusion_reaction <fix_diffusion>` is
applied, the fix also update substrate utilisation (reaction) rates R at each affected grid cell using the following
equations:

.. math::

   R_{sub} & =  \mu \cdot (\frac{1}{Y} \cdot a_{cat} + a_{ana}) \cdot X   &  \text{if } q_{cat} > \alpha \cdot m_{req}

   R_{sub}  & = \mu  \cdot a_{cat} \cdot X &  \text{if } \beta \cdot m_{req} \le q_{cat} \le \alpha \cdot m_{req}

   R_{sub}  &= - \mu \cdot a_{decay} \cdot X  & \text{if } q_{cat} < \beta \cdot m_{req}

where:

 * :math:`a_{cat}`, :math:`a_{ana}`, and :math:`a_{decay}` are the catabolic, anabolic, and decay coefficients defined in the data *file*,
 * :math:`X` is the biomass density in grid cell

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
:doc:`grid_style chemostat <grid_style_chemostat>` command.
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
*sub-i* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`.
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
*sub-i* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>` command.


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
*sub-i* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`.

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
*sub-i* is the substrate ID defined in :doc:`grid_style chemostat <grid_style_chemostat>`.

----------

.. _Gogulancea19:

**(Gogulancea, V., et al 2019)** Gogulancea, V., et al.,
Individual Based Model Links Thermodynamics, Chemical Speciation and
Environmental Conditions to Microbial Growth, Frontiers in Microbiology (2019)

.. _cabaleiro15:

**(González-Cabaleiro, R., et al 2015)** González-Cabaleiro, R., et al.,
Microbial catabolic activities are naturally selected by metabolic energy harvest rate,
ISME J (2015)