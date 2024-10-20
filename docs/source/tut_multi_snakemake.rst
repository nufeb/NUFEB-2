Managing Multiple Simulations with Snakemake
============================================

.. _Snakemake: https://snakemake.readthedocs.io/en/stable/ 

By the end of the tutorial, you will have been given an example of how
to automate and reproducibly generate and run multiple NUFEB
simulations.

An experiment will rarely rely upon a single NUFEB simulation, given
both the stochastic nature of the simulations and the number of
experimental factors likely varied.

Although it is possible to manually set up and run each NUFEB
simulation, it quickly becomes inefficient. It is also error prone and
harder to reproduce.

Various solutions exist, including writing custom scripts (*e.g.,* Bash or
Python) or the use of workflow managers. Here, we will explain how to
use Snakemake_, a popular workflow manager, to accomplish
automated, reproducible batches of NUFEB simulations.

Overview
--------

NUFEB is stochastic in nature and we would like to see how re-running
the same simulation with different random events affects the results. We
are also interested in how the initial number and spacing between
bacteria affect results.

Two major processes are involved:

1. Generating a simulation
2. Running the simulation

As part of this tutorial, you will learn how to perform those processes and to
express them as Snakemake rules.

Why Snakemake
-------------

Within each simulation, there are dependent steps. The inputfile specifying 
parameters unique to the simulation must be created before the simulation
can be run, for example.

Snakemake is an excellent tool for handling such tasks in a way which is
resource efficient, automated, and re-entrant.

By resource efficient, we mean that Snakemake can infer the order
of complex, dependent tasks. If task A is still running but task B doesn’t
require anything from task A, Snakemake knows it can start task B before
A is completed.

By automated, we mean that once a Snakemake ``Snakefile`` script is prepared, 
managing the simulations is largely hand’s off. As a bonus, the script serves 
as exact documentation of how the simulations were managed and, as a well-used
system, a Snakemake script is likely to be better understood by more users 
than a custom set of workflow scripts.

By re-entrant, imagine cancelling a set of simulations unexpectedly.
Snakemake can infer what work was completed and what work remains,
allowing you to restart the set of simulations without re-running those
which had already completed.

Essential resources
-------------------

In this case, to run multiple simulations we will need:

1. A template simulation from which others can be generated
2. Metadata concerning the variations we wish to run, which here are a
   list of random seeds and a set of grid spacings
3. Scripts which use the template and metadata to generate NUFEB
   simulation projects
4. Scripts which can initiate an individual simulation
5. A compiled, working NUFEB executable
6. A Snakemake ``Snakefile`` which coordinates the activities of all the
   above

Apart from the NUFEB executable, all of the above are included in the 
``examples/snakemake-example`` directory.

The structure within that directory is based on Snakemake recommendations, 
and the essential directories are:

1. The root directory ``./`` which will contain the ``Snakefile``
2. A
   ``./resources`` directory containing the files specifying the metadata. It 
   also contains a subdirectory, ``./resources/template`` containing the 
   template simulation
3. A ``./scripts`` directory which will contain the scripts relevant to
   generating the simulations
4. A ``./results`` directory which will contain subdirectories, each
   representing an individual simulation, its results, and the scripts
   used to initiate it

The NUFEB executable is not part of this directory tree and
should reside elsewhere. There is no need for a unique executable for
each simulation. 

.. Attention:: The location of your NUFEB executable must be edited into the ``Snakefile`` for this tutorial to work.

There is no iron-clad convention to this directory structure or the names used, 
all relevant configuration information will be in the ``Snakefile`` and can be 
altered as desired.

Setting up and running multiple simulations in Snakemake
--------------------------------------------------------

This section goes step-by-step through the example ``Snakefile`` at a high 
level of detail; we explain how rules are constructed and can call scripts,
but do not discuss the script implementations themselves. For those interested,
implementation details are discussed in a later section.


We will also attempt to explain Snakemake itself well-enough that a reader 
unfamiliar with it can understand how a rule works, but we do refer any 
confused readers to the official Snakemake tutorials.

.. seealso::
   
   `The official Snakemake tutorials <https://snakemake.readthedocs.io/en/stable/tutorial/basics.html>`_

   `The Carpentries Getting Started with Snakemake activity <https://carpentries-incubator.github.io/workflows-snakemake>`_


The basics of Snakemake
^^^^^^^^^^^^^^^^^^^^^^^

Snakemake is, at its heart, an engine which is informed by a set of
rules. Each rule specifies required inputs, expected outputs, and some
procedure which transforms those inputs into outputs. 

The rules are often abstract; conceptually there may a be a rule “say_a_name” 
and, given a list of names, Snakemake will know it will have to run the rule
once for each name. 

Given those rules, an expected endpoint, and the current state of the files, 
Snakemake can determine which rules need to be run (or re-run) and will then 
automatically execute, in the correct order, the procedures specified by each 
rule.

Specifying variables for use within snakemake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate and run simulation variations, Snakemake needs some way to
know what those variations are.

Here, those variations are the random seed to use and the arrangement of
initial cells (total number and spacing). Because these are fairly
simple variations and because Snakemake is based in Python, we can use
simple python code before defining any rules to setup our variables.

First, we create a list of random seeds, based on a text file in
``./resources`` which contains each seed on its own line.

.. highlight:: python

::

   RANDS = []                                                                      
   with open('resources/seedlist.txt', 'r') as f:                                  
       for rand_seed in f:                                                         
           RANDS.append(rand_seed.strip())                                         

Then, we read a ``json`` file which contains grid sizes and spacings, and
assign those values to variables. We also define the initial bacterial
diameter, which in this set of simulations is always 1 micron.

::

   with open('resources/grids.json','r') as f:                                     
       grids = json.load(f)                                                        
                                                                                   
   grid_N = grids['grid_N']                                                        
   grid_spacing = grids['grid_spacing']                                            
   bug_diameter = 1e-6                                                             

Finally, we indicate the location of the NUFEB executable and the number
of cores on which to run each simulation.

::

   nufebex = #insert path to NUFEB executable here, it may also be useful to access it via an ENV variable
   nufebcores = 1   

Specifying the final set of output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Snakemake works by beginning with the end in mind; it looks for the ultimate 
output, and works its way backwards to infer which rules must be run. 
The top-level rule is traditionally called ‘all’ and is used to infer the final 
desired output files. 

Internally, Snakemake identifies rules which create those final files
(if they don’t exist) and adds them to its list. It then continues iteratively,
looking at inputs for the newly found rules, finding associated rules for
those inputs, adding them to its list, and continuing on until all discovered inputs have
associated existing files or rules which will produce those files.

If a rule doesn’t contribute to the chain, it doesn’t get executed. 
Additionally, a rule doesn’t get executed if its output already exists and is
up to date.

Let’s look at our top-level rule and examine it in some detail. Many basics of
Snakemake are explained here and will be useful for understanding
subsequent rules:

::

   rule all:                                                                       
       input:                                                                      
          expand("results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/done.tkn", 
                      N = grid_N, space = grid_spacing, randno = RANDS)

The first line ``rule all:`` defines the rule name.

The next line ``input:`` specifies that we will be declaring the
files this rule will be processing. Because this is the top-level
rule, it also effectively defines the final expected outputs.

You’ll note we don’t explicitly list every unique file. Rather, we
declare the general form of the filenames, indicate the parts that vary, and
tell Snakemake how to fill in that information. 
In this case, assume we have two possible N’s (3,4), one spacing (3.5),
and two random seeds (1701, 42). Explicitly, this would be the expected
filenames:

::

   .results/3x3_3.5_default_mu_ks_yield_conc/rand1701/done.tkn, 
   .results/3x3_3.5_default_mu_ks_yield_conc/rand42/done.tkn, 
   .results/4x4_3.5_default_mu_ks_yield_conc/rand1701/done.tkn, 
   .results/4x4_3.5_default_mu_ks_yield_conc/rand42/done.tkn 

We use names in curly braces to indicate placeholders for the parts of
these names which vary, for example, if only the random seed varied, we
could call it ``{rando}``. These `wildcards <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards>`_ are a key feature of Snakemake.

The list would then be shortened to:

::

   .results/3x3_3.5_default_mu_ks_yield_conc/rand{randno}/done.tkn, 
   .results/4x4_3.5_default_mu_ks_yield_conc/rand{randno}/done.tkn 

Similarly, if we replace the square grid dimension with ``{N}`` and the
spacing with ``{space}``, we arrive at the first argument to the
``expand`` call above, namely:

::

      results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/done.tkn

What ``expand`` does is provide a mapping of the variables we defined
before to the patterns we just specified. When there is more than
one pattern, ``expand``, by default, produces all possible
combinations. There are other ways to control this, as described in the
Snakemake `explanation of expand <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function>`_.

Here, the call to ``expand`` is:

::

    expand("results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/done.tkn", 
                      N = grid_N, space = grid_spacing, randno = RANDS)

It can be read as ’Fill in ``{N}`` with values from ``grid_N``,
``{space}`` with the values of ``grid_spacing``, and ``{randno}`` with all
the values of ``RANDS``.

Specifying the rule which runs simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file required by the top-level rule, in this case, is a blank file
named ``done.tkn`` which is only generated when the simulation completes
successfully. So the next logical rule to define is the one which runs
the simulations and generates ``done.tkn`` files.

::

   rule run_sims:                                                                  
       input:                                                                      
           "results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/inputscript.nufeb"
       output:                                                                     
           "results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/done.tkn"
       run:                                                                        
           shell("results/{wildcards.N}x{wildcards.N}_{wildcards.space}_default_mu_ks_yield_conc/rand{wildcards.randno}/Allclean.sh")
           shell("results/{wildcards.N}x{wildcards.N}_{wildcards.space}_default_mu_ks_yield_conc/rand{wildcards.randno}/Allrun.sh {nufebex} {nufebcores}")

The rule, named ``run_sims`` specifies that it requires as input an
``inputscript.nufeb`` file in a directory whose name follows a pattern based on
the grid size, bacteria spacing, and random number.  Note that you do not have to use
``expand`` here, as that was already taken care of by ``rule all:`` and
those changes, conceptually, ‘trickle down’.

The rule also specifies that upon running it will generate an output, a
``done.tkn`` file in the same pattern-matched directory. Note how this
output matches the pattern of the inputs in ``rule all:``, this is core
to how Snakemake infers chains of rules.

Finally, the rule specifies two shell commands to run in order for each
specific simulation. First, a script named ``Allclean.sh`` will be
run. Its purpose is to remove any simulation run data if the directory
already exists. It is probably not necessary here, as the directory
shouldn’t exist, but it is good practice for every NUFEB simulation to
have a defined cleanup script and to run that before any simulation.

Second, a script named ``Allrun.sh`` is run. This is the script which
handles the actual call to the NUFEB executable.

Both scripts live in the pattern-matched directories listed above and
are accessed via the ``shell`` command. Because of how Snakemake works,
the varying parts of the pattern must
be prefixed here with ``wildcard``.

Notice also how ``Allrun.sh`` has two arguments, the path to the
executable and the number of cores to use. These are specified by
telling Snakemake to to use the ``nufebex`` and ``nufebcores``
variables. These do not need the ``wildcard`` prefix.

.. note:: The ``done.tkn`` file can be generated automatically in an inputscript by including this as the last line: ``shell touch done.tkn``

Specifying the rule which generates simulation specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As you might infer, the ``run_sims`` rule requires a simulation
specification given by ``inputscript.nufeb``. (NUFEB runs also often
include an ``atom.in`` file, but because both are generated at the same
time, the ``atom.in`` file was omitted for brevity.)

The following rule generates the simulation specification:

::

   rule gen_sims:                                                                  
       input:                                                                      
          "resources/seedlist.txt",                                                
          "resources/grids.json",                                                  
          "resources/template/inputscript.nufeb",                                    
          "resources/template/Allrun.sh",                                          
          "resources/template/Allclean.sh"                                         
       output:                                                                     
          "results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/inputscript.nufeb",
          "results/{N}x{N}_{space}_default_mu_ks_yield_conc/rand{randno}/atom.in"  
       run:                                                                        
          shell("python scripts/generate.py {wildcards.randno} {wildcards.N}x{wildcards.N}_{wildcards.space}_default_mu_ks_yield_conc/ {wildcards.N} {wildcards.N} {wildcards.space} 1e-6 2e-4")
          shell("python scripts/create_grid.py {wildcards.randno} {wildcards.N}x{wildcards.N}_{wildcards.space}_default_mu_ks_yield_conc/ {wildcards.N}           {wildcards.space}")

This rule declares it will generate the two output files
(``inputscript.nufeb`` and ``atom.in``) which together specify a unique NUFEB
simulation. It will do so in a now-familiar set of pattern-matched
directories.

The rule also states that it requires a number of input files.
Notably, these are files which we’ve created by hand rather than expect
to be created as part of Snakemake. This is a hallmark of an ‘entry
point’ rule. Essentially these input files are either metadata
explaining the variations or a template simulation description to base
all concrete descriptions upon.

Finally, the rule states that it will run two python scripts (again, via
``shell``). The first script, located in ``.scripts/generate.py`` is
responsible for generating an ``Inputscript.lmp`` and takes a number of
arguments, such as the target directory in which to create the file, the
number of grid rows and columns, the spacing between bacteria, the
default bacteria size (1e-6) and the height of the simulation box
(2e-4). Notice that the final two arguments are constant throughout all
simulations and can be specified as plain values.

The second script generates the associated ``atom.in`` file and takes
similar arguments.

Running it all
^^^^^^^^^^^^^^

One of the simplest invocations is running  
``snakemake -j <ncores> Snakefile`` in the top-level directory where the
``Snakefile`` resides. This tells Snakemake to infer the rule chain
based on the information in ``Snakefile`` and to execute the rules.
As long as specific instances of rules do
not depend on each other, Snakemake can automatically run them in
parallel.

You may wish to first use a dry-run with the ``-n`` flag:
``snakemake -n -j <ncores> Snakefile`` This causes to Snakemake to infer
the rule chain and tell you what it plans to do, but it does not execute
any of the rules themselves. It is very useful for quickly seeing if the
ruleset does what you intend to do.

Implementation details
----------------------

These sections describe, for those interested, implementation details of
scripts and metadata used by the Snakemake rules. As with much code,
there are many possible implementations, as such, these implementations
are meant to describe one way to go about things, not to prescribe THE
way to do things. 

Creating a simulation template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation template is simply a *nearly* valid NUFEB inputscript. In
two places, we’ve replaced exact values with keywords which are not
valid for an inputscript but which will be replaced by the simulation
generation code. There is no set list of keywords, they are just unique
identifiers which we will write a script to recognize.

In the first case we use ``het_seed`` to indicate where we want to
insert each simulation’s random seed:

::

   # biological model fixes                                                        
   fix div HET_ALL nufeb/division/coccus 1.36e-6 het_seed 

In the second case, we use ``template_fullvol`` to indicate where we
want to insert the simulation volume, which varies based on the initial number
and spacing of organisms.

::

   variable fullvol equal "template_fullvol" 

Note how we take advantage of LAMMPS’ ``variable`` keyword to reduce the
number of places in which we need to insert the template keyword.

You may have noticed the template simulation does not include the
``atom.in`` file which usually denotes specified the initial bacteria
present in the system. That is because in this case it will be entirely
generated from the metadata.

Simulation generation scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two simulation generation scripts are used in combination to generate
each complete unique simulation specification (*i.e.,* a valid NUFEB
inputscript and associated ``atom.in``). Both scripts take 
arguments which will be populated by Snakemake. The location of the
template directory is hardcoded but could have been made into another
argument if such flexibility was desired.

Generating the NUFEB inputscript
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main steps are to:

1. specify a nearly complete inputscript with placeholders indicating values which will change between simulations
2. based on input, alter the inputscript to specific values
3. save the alterations to an individual simulation directory


The NUFEB inputscript is generated by ``generate.py``, which takes 7
arguments: 

1. the specific number used to seed the random number generator 
2. the directory name (relative to ``./results``) in which to create the simulation directory

The subesquent arguments are used to calculate the simulation volume and
are: 

3. the number of bacterial rows 
4. the number of bacterial columns
5. the spacing between bacteria 
6. the initial diameter of a bacterium
7. the height of the simulation box

Obviously, some aspects of the script are unique to the experiment. The
essential part consists of: 

1. reading the template inputscript 
2. replacing any instances of ``het_vol`` and ``template_fullvol`` with values of the first argument and the calculated simulation volume (``box_vol``) 
3. saving the results as an inputscript in the correct simulation directory (``targetdir``)

We can accomplish this with ``fileinput`` as below:

::

   fname = os.path.join(targetdir, 'Inputscript.lmp')                              
       with fileinput.FileInput(fname, inplace=True, backup='.bak') as file:           
           for line in file:                                                           
               het_line = line.replace('het_seed', sys.argv[1])                        
               print(het_line.replace('template_fullvol', str(box_vol)), end='') 

For more intensive token replacement, you may want to look into a
purpose built templating engine, such as `Jinja <https://jinja.palletsprojects.com/en/3.1.x/>`_. Also, Python
is not the only option. 

Generating the initial organisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The layout and characteristics of the initial microbial population is
specified in a separate file for NUFEB simulations, often called
``atom.in``.

Here, we generate the ``atom.in`` using ``create_grid.py``, which takes
4 arguments. 

1. the specific number used to seed the random number generator 
2. the directory name (relative to ``./results``) in which to create the simulation directory 
3. the number of bacterial rows/columns (assumed square grid) 
4. the horizontal/vertical distance between each bacteria

From this, the script generates an atom file describing an NxN grid of
bacterial cells with the specified spacing. All other aspects of the
cells are given default values.

Further reading
---------------

.. seealso::
    `The paper which necessitated this approach <https://royalsocietypublishing.org/doi/full/10.1098/rsfs.2023.0010>`_
        Weaver, J. E. (2023). Quantifying drift-selection balance using an agent-based biofilm model of identical heterotrophs under low-nutrient conditions. Interface Focus, 13(4), 20230010. doi:`https://doi.org/10.1098/rsfs.2023.0010 <https://doi.org/10.1098/rsfs.2023.0010>`_
        
    `Offical Snakemake Documentation <https://snakemake.readthedocs.io/en/stable/index.html>`_
        I've glossed over many details in Snakemake and played a little loose with terminology. There is a *a lot* that Snakemake can do for you and it's worth learning.
    
    `Example directory <https://github.com/nufeb/NUFEB-2/tree/master/examples/snakemake-example>`_
        Location of the relevant example in the repository. 
