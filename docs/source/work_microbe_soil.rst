Workshop: modelling microbial dynamics in soil  
==================================================

This workshop aims to use NUFEB to simulate the growth and mobility of microbes in soil environment.

Run simulations on Code Ocean
--------------------------------

1. Create a Code Ocean account at https://codeocean.com/
2. Email your account (Email address) to bowen.li2@newcastle.ac.uk so that we can share the code to you.
3. You should be able to see a new Capsule named NUFEB in your profile. Enter the Capsule and create a copy of it by clicking the menu **Capsule -> Duplicate**.
4. The NUFEB capsule provides two examples in the ``/data`` directory for modelling microbes in soil. 
   ``/data/soil-dry`` simulates microbial growth on soil particle surface in dry condition, and ``/data/soil-sat`` simulates microbial growth and mobility in water saturated condition.
5. To run one of them, go to the ``run-example`` file and choose the example that you want to run by uncomment corresponding line. 
   Then right click on the file in the left *File* panel, and select **Set as File to Run**.
6. Click the **Reproducible Run** button in the upper-right corner. During the simulation, information such as total number of individuals and biomass over time are printed in the console panel, 
   and output files are saved in the ``/results`` directory.
7. To download output files, right click on the ``/results`` directory and select **Download**.


Run NUFEB on Code Ocean
--------------------------------

In this section, we will extend the ``/data/soil-dry`` case to simulate competition between two microbial types - yield strategist (YS) and rate strategist (RS). 
The former refers to a group of microbes with high growth yield at low growth rate 
which is a case of an altruistic strategy because they increase the fitness of the group 
by using resources economically at the cost of decreased fitness, or growth rate, of the individual. 
On the other hand, RS is a group of microbes with high growth rate at low growth yield which is a case of an egoistic strategy.

It was considered that RS will outcompete YS in high substrate 
concentration condition, as the high concentration can penetrate into the microbial cluster deeply (due to diffusion) that makes the overall growth rate of RS faster than YS.
While substrate-limit condition could in favour the growth of YS. This is because low concentration area is formed inside clusters, 
the economical use of common limited resources for YS increase their fitness to such environment. 

To model and simulate the YS/RS system, we assume the existing microbial type in the ``/data/soil-dry`` case is YS 
(type 2, growth rate = 2e-4, Yield = 0.43, and Ks = 2e-4). 
We need to introduce the other microbial type RS to the model 
with growth rate = 4e-4, Yield = 0.21, Ks = 4e-4. This requires to modify or add some commands in the input script: 

1. Open the file ``data/soil-dry-compete/Inputscript.lmp``, which is a copy of ``/data/soil-dry/Inputscript.lmp`` with mono-microbial type.
2. To introduce the new type RS, we need to change/add the following commands in the file:

 a. Line 28, change: 
   
  .. parsed-literal::
  
    create_box 3 box     # use 3 atom types in the box
   
 b. Line 55, add: 
 
  .. parsed-literal::
  
    # randomly place 100 type 3 initial microbes on the surface
    create_atoms 3 random 100 2345 soil_reg var v_ibac set x xx set y yy set z zz
    
 c. Line 62-63, add: 
 
  .. parsed-literal::
  
    # set attributes for type 3
    set type 3 diameter 1.5e-6       
    set type 3 density 30 
   
 d. Line 66, change:
 
   .. parsed-literal::
  
    group  bac type 2  3  #assign type 2 and 3 atoms to bac group
    
 d. Line 67-68, add:
 
   .. parsed-literal::
  
    group  YS  type 2     #assign type 2 atoms to YS group     
    group  RS  type 3     #assign type 3 atoms to RS group
    
 e. Line 96, change:
 
   .. parsed-literal::
  
    # change group-ID to YS
    fix f_grow1 YS nufeb/growth/monod sub 2e-04 growth 0.0002 yield 0.43 
    
 e. Line 97, add:
 
   .. parsed-literal::
  
    # impose microbial growth process and define kinetics parameters for RS
    fix f_grow2 RS nufeb/growth/monod sub 4e-04 growth 0.0004 yield 0.21 
    
 f. Line 117-118, add:
   
   .. parsed-literal::
   
    # assigns numbers of YS and RS to corresponding variables 
    variable nYS equal "count(YS)"
    variable nRS equal "count(RS)"
   
 g. Line 120, change:
 
    .. parsed-literal::
    
     # output YS and RS population sizes
     thermo_style custom step cpu v_nbac v_mass v_nYS v_nRS
 
3. Simulate the model by uncomment the line *case=soil-dry-compete* in the ``run-example`` file. Check the numbers of RS and YS at the end of the run.

4. In ``data/soil-dry-compete/Inputscript.lmp``, change the initial concentration from 1e-4 (kg/m3) to 1e-5 by altering grid_modify (Line 74) command to:

    .. parsed-literal::
    
     grid_modify  set sub  pp pp dn 0 1e-5
    
   Also change the total simulation time to 220 (220 hours)
    
    .. parsed-literal::

     run 220
   
   Rerun the simulation and see the difference of RS and YS populations at the end of the run comparing to the previous result.
 