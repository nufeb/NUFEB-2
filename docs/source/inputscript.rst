Input Script
==============

In order to run NUFEB simulation, an inputscript (text file) is usually prepared with
certain commands and parameters list. NUFEB executes by reading those commands and
parameters, one line at a time. When the input script ends, NUFEB exits. Each command
causes NUFEB to take some actions. It may set an internal variable, read in a file, or run
a simulation. Many input script errors are detected by NUFEB and an ERROR or WARNING message is printed. 
The documentation for each command lists restrictions on how the command can be used.

Input script structure
------------------------

This section describes the structure of a typical NUFEB input script. We will take the
script in NUFEB-dev/examples/biofilm-het directory as an example for the explanation.

A NUFEB input script typically consists of the following parts:

1. :ref:`General settings <settings>`
2. :ref:`Microbes, nutrients and simulation box <init>`
3. :ref:`IbM processes <ibm>`
4. :ref:`Computations and outputs <output>`
5. :ref:`Run a simulation <run>`


.. _settings:

General settings
------------------------------

Set parameters that need to be defined before atoms and substrates are created or read-in from a file.


The relevant commands are 
:doc:`atom_style coccus <atom_vec_coccus>`,
:doc:`atom_style bacillus <atom_vec_bacillus>`,
`units* <https://docs.lammps.org/units.html>`_,
`newton* <https://docs.lammps.org/newton.html>`_,
`processors* <https://docs.lammps.org/processors.html>`_, 
`atom_modify* <https://docs.lammps.org/atom_modify.html>`_,
`comm_modify* <https://docs.lammps.org/comm_modify.html>`_.
(Commands with asterisk `*` are from LAMMPS without modifications.)

Example:

.. code-block:: 

	#-------------------------------------General settings--------------------------------------#
	
	units si                                    # using si units
	
	atom_style      coccus                      # using nufeb coccus atom style
	atom_modify     map array sort 1000 1e-6    # map array: find atoms using indices
                                                    # sort 1000 1e-6: sort every 
                                                    # 1000 steps with 1e-6 binsize
	
	newton          off                         # forces between local and ghost atoms are 
	                                            # computed in each processor without communication
	
	processors      * * 1                       # processor grid
	
	comm_modify     vel yes                     # communicate velocities for ghost atoms
	comm_modify     cutoff 2e-6                 # guarantee that enough atoms are
	                                            # communicated to correctly compute




.. _init:

Microbes, nutrients and simulation box 
---------------------------------------

Define initial microbes, nutrients (substrates) that can be metabolised by the microbes, 
and simulation box that microbes reside 

Substrates are defined via grid_style command 
(:doc:`grid_style chemostat <grid_style_chemostat>`).
There are 3 ways to define simulation box and initial microbes in NUFEB. 
Read them in from
(1) a data file via the `read_data* <https://docs.lammps.org/read_data.html>`_ command,
(2) a restart file via the 
`read_restart* <https://docs.lammps.org/read_restart.html>`_ commands,
or (3) create a simulation cell and fill it with atoms on
a lattice, using these commands:
`lattice* <https://docs.lammps.org/lattice.html>`_, 
`region* <https://docs.lammps.org/region.html>`_, 
`create_box* <https://docs.lammps.org/create_box.html>`_, 
`create_atoms* <https://docs.lammps.org/create_atoms.html>`_ or
`read_dump* <https://docs.lammps.org/read_dump.html>`_.

Example:

.. code-block:: 

	#---------------------Microbes, nutrients and simulation box----------------------------#
	
	boundary        pp pp ff                   # periodic domain boundaries in x and y
	                                           # fixed boundary in z
	                                            
	read_data       atom.in	                   # read atom.in file which defines box size
	                                           # and initial atoms	                
	
	group           het   type 1               # assign type 1 atoms to het group
	group           eps   type 2               # assign type 2 atoms to eps group
	
	neighbor        5e-7 bin                   # setting neighbour skin distance and style
	neigh_modify    check yes                  # rebuild neighbour list if any atom
	                                           # had moved more than half the skin distance
	
	# use nufeb/chemostat grid style, define substrate types and diffusion grid size
	grid_style      nufeb/chemostat 4 sub o2 no2 no3 4e-6  
	
	# set diffusion boundary conditions and initial concentration 
	grid_modify     set sub  pp pp nd  1e-4 1e-4
	grid_modify     set o2   pp pp nd  1e-4 1e-4
	grid_modify     set no2  pp pp nd  1e-4 1e-4
	grid_modify     set no3  pp pp nd  1e-4 1e-4


.. _ibm:

IbM processes
------------------------------

NUFEB provides a variety of individual-based microbial modelling (IbM) processes.
They are classified into different submodules depending on the timesteps and 
their execution orders in the NUFEB integration procedure 
(:doc:`run_style nufeb <run_style_nufeb>` command):

1. :doc:`Biological processes <list_biology>` 
2. :doc:`Physical processes <list_physics>` 
3. :doc:`Chemical processes <list_post_physics>`
4. :doc:`Chemical processes <list_chemistry>`
5. :doc:`Reactor processes <list_reactor>`

Example:

.. code-block:: 

	#-----------------------------------Biological processes-------------------------------------#
	
	# heterotrophs growth
	fix growth_het het nufeb/growth/het sub 3.5e-5 o2 0 no2 0 no3 0 growth 0.00028 yield 0.61 decay 0.0 epsyield 0.18 anoxic 0.0 epsdens 30
	
	# heterotrophs division 
	fix div het nufeb/division/coccus 1.36e-6 30 12345
	
	# EPS extraction from heterotrophs
	fix eps_ext het nufeb/eps_extract 2 eps 1.3 30 12345
	
	
	#------------------------------------Physical processes--------------------------------------#
	
	# contact force between atoms
	pair_style  gran/hooke/history 1e-4 NULL 1e-5 NULL 0.0 1
	pair_coeff  * *
	
	# contact force between atoms and domain boundary
	fix wall all wall/gran hooke/history 1e-3 NULL 1e-4 NULL 0 0 zplane 0.0 8e-5
	
	# viscous damping force between atoms
	fix vis all viscous 1e-5
	
	# NVE integration with maximum distance limit
	fix nve all nve/limit 1e-8
	
	#-----------------------------------Post-physical processes----------------------------------#
	
	# dynamic diffusion coefficient based on biofilm density
	fix coeff_sub all nufeb/diffusion_coeff sub ratio 0.8
	
	#-----------------------------------Chemical processes---------------------------------------#
	
	# diffusion reaction for substrate
	fix diff_sub all nufeb/diffusion_reaction sub 1.6e-9 
	
	
.. _output:

Computations and outputs
------------------------------

Define :doc:`computations <list_compute>`  to compute various microbial properties during a simulation.

NUFEB supports output simulation date in the following formats: 

1. Plain text (`thermo_style* <https://docs.lammps.org/thermo_style.html>`_ and  `thermo* <https://docs.lammps.org/thermo.html>`_ commands)
2. VTK (:doc:`dump vtk/grid <dump_vtk_grid>` and `dump vtk* <https://docs.lammps.org/dump_vtk.html>`_ commands)
3. HDF (:doc:`dump hdf5 <dump_hdf5>` command)
4. PNG or JPEG or PPM images, or as a single movie file (`dump_image* <https://docs.lammps.org/dump_image.html>`_ command)

Example:

.. code-block:: 

	#----------------------------------Computations and Outputs---------------------------------#
	
	# compute biofilm pressure
	compute vol all nufeb/volume
	compute ke all ke
	variable one equal 1.0
	compute press all pressure NULL pair vol v_one
	variable press equal "(c_ke + c_press) / (3.0 * c_vol)" 
	
	# compute total mass
	variable mass equal "mass(all)"
	
	# dump atom and grid data to /vtk folder in vtk format
	shell mkdir vtk
	dump du1 all vtk 10 vtk/dump*.vtu id type diameter
	dump du2 all grid/vtk 10 vtk/dump_%_*.vti con rea den gro
	
	# screen and log outputs
	thermo_style custom step cpu atoms v_press v_mass
	thermo 1

.. _run:

Run
------------------------------

A NUFEB simulation is run using the :doc:`run_style nufeb <run_style_nufeb>` and `run* <https://docs.lammps.org/run.html>`_ commands.

Example:

.. code-block:: 

	#--------------------------------------Run------------------------------------------------#
	
	# issue run command, define timesteps for physical and chemical processes
	run_style nufeb diffdt 1e-4 difftol 1e-6 pairdt 1e-2 pairtol 1 pairmax 1000 diffmax 5000
	
	# define biological timesteps
	timestep 1000
	
	# run 900 biological steps (9x10^5 seconds)
	run 900
