Commands
============

This section lists all NUFEB commands that can be used
for (individual-based) modelling of microbial systems. 
Note that commands with asterisk `*` are inherited from LAMMPS with/without modifications.
Some commands are part of NUFEB optional packages,
which means they cannot be used unless the packages 
was included when building NUFEB. 

.. contents:: 
		:local:
		:depth: 1
   




.. _cmd_1:

.. _comm:


General settings
-------------------------------------------

+--------------------------------------------+---------------------------------------+
| :doc:`grid_style monod <grid_style_monod>`: grid style for Monod-based growth      |
+--------------------------------------------+---------------------------------------+
| :doc:`grid_style reactor <grid_style_reactor>`: grid style for reactor coupling    |
+--------------------------------------------+---------------------------------------+

Microbe
-------------------------------------------

+--------------------------------------------+-----------------------------------------------+
| :doc:`atom_style coccus <atom_vec_coccus>`: spherical microbe                              |
+--------------------------------------------+-----------------------------------------------+
| :doc:`atom_style bacillus <atom_vec_bacillus>`: rod-shaped microbe                         |
+--------------------------------------------+-----------------------------------------------+
| :doc:`fix nufeb/property/cycletime <fix_property_cycletime>`: cell cycle time property     | 
+--------------------------------------------+-----------------------------------------------+
| :doc:`fix nufeb/property/generation <fix_property_generation>`: cell generation property   |
+--------------------------------------------+-----------------------------------------------+
| :doc:`fix nufeb/property/plasmid <fix_property_plasmid>`: plasmid property                 |
+--------------------------------------------+-----------------------------------------------+


Biological processes
-------------------------------------------

+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/simple <fix_monod_simple>`: Monod-based growth model                     |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/anammox <fix_monod_anammox>`: Monod-based growth for ANAMMOX             |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/aob <fix_monod_aob>`: Monod-based growth for ammonia-oxidizing bacteria  |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/cyano <fix_monod_cyano>`: Monod-based growth for cyanobacteria           |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/ecoli/wild <fix_monod_ecoli_wild>`: Monod-based growth for ecoli (WT)    |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/eps <fix_monod_eps>`: Monod-based growth for EPS                         |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/het <fix_monod_het>`: Monod-based growth for heterotroph                 |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/monod/nob <fix_monod_nob>`: Monod-based growth for nitrite-oxidizing bacteria  |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/divide/coccus <fix_divide_coccus>`: microbe division for coccus                |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/divide/bacillus <fix_divide_bacillus>`: microbe division for bacillus          |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/divide/bacillus/minicell <fix_divide_minicell>`: abnormal division for bacillus|
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/eps_excretion <fix_eps_excretion>`: EPS excretion from heterotroph             |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/death/diameter <fix_death_diameter>`: microbe death (critical diameter)        |
+--------------------------------------------+---------------------------------------------------+
| :doc:`fix nufeb/death/plasmid <fix_death_plasmid>`: microbe death (critical # of plasmids)     |
+--------------------------------------------+---------------------------------------------------+


Physical processes
-------------------------------------------

+--------------------------------------------+------------------------------------------------------+
| :doc:`pair_style bacillus <pair_bacillus>`: pairwise interaction for bacillus                     |
+--------------------------------------------+------------------------------------------------------+
| :doc:`pair_style gran/hooke <pair_gran_hooke>`: pairwise interaction for coccus                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/adhesion <fix_adhesion>`: microbe-microbe adhesion force                          |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/adhesion/eps <fix_eps_adhesion>`: EPS-microbe adhesion force                      |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nve/bacillus/limit <fix_nve_bacillus_limit>`: constant NVE update for bacillus          |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/shear <fix_shear>`: shear force                                                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix nufeb/wall_adhesion <fix_wall_adhesion>`: wall-microbe adhesion force                   |
+--------------------------------------------+------------------------------------------------------+
| :doc:`fix wall/gran <fix_wall_gran>`: wall-microbe frictional force                               |
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

Compute  
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
| :doc:`compute nufeb/plasmid/ave_copy <compute_ave_plasmid>`: average plasmid copy number                     |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/plasmid/ave_nbirth <compute_plasmid_nbirth>`: average plasmid copy number at cell birth  |
+--------------------------------------------+-----------------------------------------------------------------+
| :doc:`compute nufeb/plasmid/copy <compute_plasmid_copy>`: plasmid copy number                                |
+--------------------------------------------+-----------------------------------------------------------------+


Post-processing  
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



Other LAMMPS commands
-------------------------------------------
This section lists all LAMMPS commands that are tested to be compatible 
with NUFEB, but they are not directly related to microbial modelling. 
