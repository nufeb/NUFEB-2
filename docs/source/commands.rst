Commands
============

This section lists all NUFEB commands and 
essential LAMMPS commands (with asterisk `*`) that can be used
for (individual-based) modelling microbial systems. 
Some commands are part of NUFEB or LAMMPS optional packages,
which means they cannot be used unless the packages 
were included when building NUFEB. 

 .. note::
 
	The section only lists the LAMMPS commands which are commonly used for microbial modelling.
	Many of the commands not listed here are also compatible with NUFEB package and their descriptions
	are available in `LAMMPS user manual <https://docs.lammps.org/Manual.html>`_.


.. contents:: 
		:local:
		:depth: 1
   




.. _cmd_1:

.. _comm:


General settings
-------------------------------------------

+--------------------------------------------+-------------------------------------------------------------+
| `comm_modify* <https://docs.lammps.org/comm_modify.html>`_: settings for inter-processor communication   |
+--------------------------------------------+-------------------------------------------------------------+
| `processors* <https://docs.lammps.org/processors.html>`_: parallel decomposition scheme                  |
+--------------------------------------------+-------------------------------------------------------------+
| `units* <https://docs.lammps.org/units.html>`_:  style of units (NUFEB supports SI only)                 |
+--------------------------------------------+-------------------------------------------------------------+

Microbes, substrates and simulation box 
-------------------------------------------

+--------------------------------------------+---------------------------------------------------------+
| :doc:`atom_style coccus <atom_vec_coccus>`: spherical microbe                                        |
+--------------------------------------------+---------------------------------------------------------+
| :doc:`atom_style bacillus <atom_vec_bacillus>`: rod-shaped microbe                                   |
+--------------------------------------------+---------------------------------------------------------+
| :doc:`fix nufeb/property/cycletime <fix_property_cycletime>`: microbe attribute: cell cycle time     | 
+--------------------------------------------+---------------------------------------------------------+
| :doc:`fix nufeb/property/generation <fix_property_generation>`: microbe attribute: cell generation   |
+--------------------------------------------+---------------------------------------------------------+
| :doc:`fix nufeb/property/plasmid <fix_property_plasmid>`: microbe attribute: plasmid                 |
+--------------------------------------------+---------------------------------------------------------+
| :doc:`grid_modify <grid_modify>`: modify grid attributes                                             |
+----------------------------------------------------+-------------------------------------------------+
| :doc:`grid_style chemostat <grid_style_chemostat>`: grid style for chemostat coupling                |
+----------------------------------------------------+-------------------------------------------------+
| :doc:`grid_style simple <grid_style_simple>`: simple grid style                                      |
+--------------------------------------------+---------------------------------------------------------+
| `atom_modify* <https://docs.lammps.org/atom_modify.html>`_: modify (system) attributes of atoms      |
+----------------------------------------------------+-------------------------------------------------+
| `boundary* <https://docs.lammps.org/boundary.html>`_: boundaries for simulation box                  | 
+----------------------------------------------------+-------------------------------------------------+
| `create_atoms* <https://docs.lammps.org/create_atoms.html>`_: create atoms on a lattice              | 
+----------------------------------------------------+-------------------------------------------------+
| `create_box* <https://docs.lammps.org/create_box.html>`_: create a simulation box based on region    | 
+----------------------------------------------------+-------------------------------------------------+
| `group* <https://docs.lammps.org/group.html>`_: group a collection of atoms                          | 
+--------------------------------------------+---------------------------------------------------------+
| `lattice* <https://docs.lammps.org/lattice.html>`_: define a lattice stucture over simulation box    | 
+--------------------------------------------+---------------------------------------------------------+
| `read_data* <https://docs.lammps.org/read_data.html>`_: read external data file                      |
+----------------------------------------------------+-------------------------------------------------+
| `region* <https://docs.lammps.org/region.html>`_: defines a geometric region of space                |
+----------------------------------------------------+-------------------------------------------------+
| `neighbor* <https://docs.lammps.org/neighbor.html>`_: parameters for neighbor lists                  |
+----------------------------------------------------+-------------------------------------------------+
| `neigh_modify* <https://docs.lammps.org/neigh_modify.html>`_: modify parameters for neighbor lists   |
+----------------------------------------------------+-------------------------------------------------+

Biological processes
-------------------------------------------

+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/anammox <fix_growth_anammox>`: growth for anammox bacteria                  |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/aob <fix_growth_aob>`: growth for ammonia-oxidizing bacteria                |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/cyano <fix_growth_cyano>`: growth for cyanobacteria                         |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/ecoli <fix_growth_ecoli>`: growth for E.coli                                |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/eps <fix_growth_eps>`: extracellular polymeric substances decay             |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/het <fix_growth_het>`: growth for heterotroph                               |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/monod <fix_growth_monod>`: monod-based growth                               |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/nob <fix_growth_nob>`: growth for nitrite-oxidizing bacteria                |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/growth/simple <fix_growth_simple>`: simple linear growth                           |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/division/coccus <fix_divide_coccus>`: microbe division for coccus                  |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/division/bacillus <fix_divide_bacillus>`: microbe division for bacillus            |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/division/bacillus/minicell <fix_divide_minicell>`: abnormal division for bacillus  |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/eps_extract <fix_eps_extract>`: EPS excretion from heterotroph                     |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/death/diameter <fix_death_diameter>`: microbe death (critical diameter)            |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/plasmid/conjugate <fix_plasmid_conjugate>`: conjugation for gram-positive bacteria |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/plasmid/kill <fix_plasmid_kill>`: microbe death (critical # of plasmids)           |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/plasmid/partition <fix_plasmid_partition>`: type II partitioning system            |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/plasmid/replicate <fix_plasmid_replicate>`: mass-based plasmid replication         |
+--------------------------------------------+-------------------------------------------------------+


Physical processes
-------------------------------------------

+--------------------------------------------+------------------------------------------------------+
| :doc:`pair_style bacillus <pair_bacillus>`: pairwise interaction for bacillus                     |
+--------------------------------------------+------------------------------------------------------+
| :doc:`pair_style gran/hooke <pair_gran_hooke>`: pairwise interaction for coccus                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/adhesion <fix_adhesion>`:  adhesion force                                         |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/adhesion/eps <fix_adhesion_eps>`: EPS adhesion force                              |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/adhesion/bacillus <fix_adhesion_bacillus>`: adhesion force for bacillus           |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nve/bacillus/limit <fix_nve_bacillus_limit>`: constant NVE update for bacillus          |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/shear <fix_shear>`: shear force                                                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/wall_adhesion <fix_wall_adhesion>`: wall-microbe adhesion force                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix wall/gran <fix_wall_gran>`: wall-microbe frictional force                               |
+--------------------------------------------+------------------------------------------------------+

Post-physical processes
-------------------------------------------

+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/diffusion_coeff <fix_diffusion_coeff>`: dynamic diffusion coefficient             |
+--------------------------------------------+------------------------------------------------------+


Chemical processes
-------------------------------------------

+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/diffusion_reaction <fix_diffusion>`: substrate diffusion and reaction              |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/gas_liquid <fix_gas_liquid>`: gas liquid transfer                                  |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/reactor/gas_balance <fix_reactor_gas_balance>`: mass balance in gas phase          |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/reactor/solute_balance <fix_reactor_solute_balance>`: mass balance in solute phase |
+--------------------------------------------+-------------------------------------------------------+


Post-chemical processes
-------------------------------------------

+--------------------------------------------+-------------------------------------------------------+
| :doc:`fix nufeb/boundary_layer <fix_boundary_layer>`: dynamic diffusion boundary layer             |
+--------------------------------------------+-------------------------------------------------------+


Computes
-------------------------------------------

+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/ave_conc <compute_ave_conc>`: average substrate concentration                            |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/ave_length <compute_ave_length>`: average microbe length                                 |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/density <compute_density>`: biomass density                                              |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/volume <compute_volume>`: total microbe volume                                           |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/plasmid/ave_copy <compute_ave_copy>`: average plasmid copy number                        |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/plasmid/ave_nbirth <compute_plasmid_nbirth>`: average plasmid copy number at cell birth  |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/plasmid/copy <compute_plasmid_copy>`: plasmid copy number                                |
+--------------------------------------------+-----------------------------------------------------------------+


Outputs
-------------------------------------------

+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump image <dump_image>`: dump JPEG, PNG or PPM image files                                  |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump movie <dump_movie>`: dump movie file                                                    |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump modify <dump_modify>`: modify parameters of dump command                                |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump vtk <dump_vtk>`: dump microbe data in VTK format                                        |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump vtk/grid <dump_vtk_grid>`: dump grid data in VTK format                                 |
+--------------------------------------------+-------------------------------------------------------+
| :doc:`dump hdf5 <dump_hdf5>`: dump data in hdf5 format                                             |
+--------------------------------------------+-------------------------------------------------------+


Run
-------------------------------------------

+----------------------------------------------------+---------------------------------------+
| :doc:`run_style nufeb <run_style_nufeb>`: time integrator for NUFEB simulation             |
+----------------------------------------------------+---------------------------------------+

