Stopping simulations at a consistent  percent biomass
=====================================================

Problem
-------

You want to end simulations when the biomass occupies a specific percentage of the simulation volume, enabling consistent comparisons.

Solution
--------

Use the NUFEB :doc:`compute_volume` to calculate the biomass percentage, which the LAMMPS `halt fix <https://docs.lammps.org/fix_halt.html>`_ uses to determine if the simulation should stop.

::

    # Compute the volume of all bacteria (could be adapted for specific groups)
    compute biomass_vol all nufeb/volume

    # Convert to a percent of simulation volume 
    #(sim volume is either explicitly entered or generated via script)
    variable sim_vol equal "1e-13"
    variable biomass_pct equal "c_biomass_vol/v_sim_vol"

    #... skip to just before run_style

    # halt the simulation when the biomass reaches 10 percent
    fix halt_vol all halt 1 v_biomass_pct > 0.10 error soft

Discussion
----------

Biofilms may grow at very different rates depending on simulation conditions; sometimes it may be inappropriate to run all simulations to a common timepoint and then compare. Rather, a more 'apples to apples' comparison may involve growing biomass to a set percentage of the simulation volume.

NUFEB has a :doc:`compute_volume` which calculates the volume of all the organisms matching the relevant group identifier. Meanwhile, LAMMPS has a built-in ``halt`` fix which can perform logical comparisons to determine when (and how) to terminate a simulation.

Generally, the volume computation should occur within the input script after all relevant growth fixes.  The halt fix should be near the end of the input script, generally just prior to the ``run_style`` command.  

.. seealso::
    Computing biomass volume
        :doc:`compute_volume`

    Conditionally halting a LAMMPS script
        `LAMMPS halt fix <https://docs.lammps.org/fix_halt.html>`_

    Creating, using, and referring to variables in LAMMPS scripts
        `variable command <https://docs.lammps.org/variable.html>`_
