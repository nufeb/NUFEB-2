.. index:: fix nufeb/T6SS/contact

fix nufeb/T6SS/contact command
==============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/T6SS/contact seed & 
            num-attackers attacker1-type attacker1-effectorID attacker1-range atacker1-cooldown &
                          attacker2-type attacker2-effectorID attacker2-range atacker2-cooldown &
                          ...
                          attacker_n-type attacker_n-effectorID attacker_n-range atacker_n-cooldown &
            num-vulns     vuln1-type vuln1-toeffector vuln1-prob intox1_group intox1_type &
                          vuln2-type vuln2-toeffector vuln2-prob intox2_group intox2_type &
                          ...
                          vuln_n-type vuln_n-toeffector vuln_n-prob intox2_group intox2_type

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* seed = number used to seed random number generator
* num-attackers = number of attackers for which a T6SS will be defined
    * attacker *n*-type = type code associated with attacker *n*
    * attacker *n*-type = ID of effector associated with attacker *n* 
    * attacker *n*-range = length of the T6SS harpoon for attacker *n*
    * attacker *n*-cooldown = cooldown between firing for attacker *n* (placeholder)
* num-vulns = number of organism types for which there is a T6SS vulnerability
    * vuln *n*-type = type code associated with vulnerable organism *n*
    * vuln *n*-toeffector = ID of the effector to which the organism is vulnerable
    * vuln *n*-prob = probability [0-1] of a hit by the effector resulting in intoxication (placeholder)
    * vuln *n*-intox_group = name of the group to place the organism into if intoxicated
    * vuln *n*-intox_type = type code of the group to place the organism into if intoxicated

Examples
""""""""

.. code-block::

   # define groups for the attacker organism, an immune organism, a vulnerable one, 
   # an attacker which uses a different effector protein, and the intoxicated state of the vulnerable organism
   group attacker type 1
   group immune type 2
   group vuln type 3
   group vuln_intoxicated typed 4
   group attacker_other type 5

   # assume that the other usual fixes for growth, death, division are set up here
   
   # for intoxication to have the expected effect, define a lysis fix applicable to the intoxicated organism
   fix lysis_vuln vuln_intoxicated nufeb/T6SS/lysis sub 2e-3 0.2

   # Define T6SS interactions, we use the & line continuation to keep things readable by humans
   fix apply_t6ss all nufeb/T6SS/contact 1701 &
       2 1 1 1.3e-6 100 &
         5 5 1.3e-6 100 &
       1 3 1 1 vuln_intoxicated 4
  
Description
"""""""""""

A Type VI Secretion System (T6SS) interaction involves two organisms, the *attacker* and the *vulnerable* organism. Attacks occur when an attacker is sufficiently close to a vulnerable organism. That distance is determined by the length of the T6SS harpoon.  When an attack occurs, the vulnerable organism is exposed to a toxic effector protein. The identify of the effector protein is what actually determines vulnerability or immunity. 

Consequentially, the T6SS fix requires defining two lists: the bacteria which contain T6SS systems and the bacteria which are vulnerable to one of the associated effector proteins.

In the above example, we have used the LAMMPS line continuation ``&`` functionality to break the fix statement into logical blocks.

The first block is one line and contains much of the 'boilerplate' of setting up a LAMMPS style fix:

.. code-block::

   fix apply_t6ss all nufeb/T6SS/contact 1701 &

Here, we specify there is a ``fix`` named ``apply_t6ss`` which applies to ``all`` atoms, is defined internally as ``nufeb/T6SS/contact`` and which uses the value of ``1701`` to initialize its random number generator.

The second block contains (here) two lines and defines the qualities related to the two T6SS-capable bacteria. In practice, there will be as many lines as there are T6SS-capable bacteria.

.. code-block::

       2 1 1 1.3e-6 100 &
         5 5 1.3e-6 100 &
 
The first number in the first line indicates that there are two T6SS-capable bacteria.  The rest of the line defines the first bacterium: the organism identified with the type code ``1`` will be associated with an effector given the id of ``1``, will have a harpoon length of ``1.3e-6`` microns, and a cooldown of 100s. There is no reason that the effector ID needs to be the same as the attacker type, but it helps to keep logical organization. The harpoon length was chosen based on the representative starting diameter of the organism type, which is in line with T6SS biology.  The cooldown parameter is currently unused by the implementation, as traditional NUFEB biological timesteps are much longer than the effective cooldown time.

The second line above defines the second bacterium, which has a type code and effector ID of ``5``

The final block defines those organisms which are vulnerable to any of the now-defined T6SS effectors and, like the previous block, can be multiple lines.
Here, it is one line:
       
.. code-block::

       1 3 1 1 vuln_intoxicated 4

The first number in the first line indicates that only one organism in this simulation is vulnerable to T6SS attacks.  That organism is the one associated with an atom type of ``3`` and it is specifically vulnerable to the effector with the ID of ``1``.  The next ``1`` indicates that there is 100% probability that an attack will result in intoxication. On intoxication, the individual will have its group changed to ``vuln_intoxicated`` and its atom type changed to ``4``.  Note that the attack success probablity is a placeholder and has no current effect on the simulation, as many T6SS attacks, particularly those studied with NUFEB, are very effective with essentially 100% probability of success especially when considering the multiple inevitable hits when staying relatively stationary within a biofilm.


The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.

The value of the ``sub-ID`` must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.

This fix requires that NUFEB be installed/built by including the ``--enable-t6ss`` argument when installing via ``install.sh`` or by including ``make yes-t6ss`` as part of the build process. 

Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`

Two parameters, the attacker cooldown and vulnerable probablity are placeholders for potential future functionality but do not affect the simulation. 

Currently, each attacker can have 1 effector defined. Bacteria can be vulnerable, however, to more than one effector.
