Computing species relative abundances
=====================================


Problem
-------

You wish to calculate at each time step the relative abundances of each species in the simulation.

Solution
--------

Use a combination of ``variable`` and ``count`` in the inputscript to compute the relevant values at each time step. The variables will then be available for subsequent inputscript commands relating to simulation control, further calculations, or output either to the screen or a file.

Assume that the relevant bacterial species are identified by groups named ``HET_1`` and ``HET_2``.

.. code-block::

    # Get the per-species counts
    variable het1_num equal "count(HET_1)"
    variable het2_num equal "count(HET_2)"

    # Calculate the combined number of both species (see note in discussion)
    variable het_all_num equal "v_het1_num + v_het2_num"

    # Calculate the relative abundances of each species
    variable het1_rel_abs equal "het1_num/het_all_num"
    variable het2_rel_abs equal "het2_num/het_all_num"

Discussion
----------

During execution of the inputscript, many values can be computed for output. The ``count`` function makes it possible to keep updated tallies.

.. tip::
   The relative abundances could also be done on either mass or volume basis using the ``mass`` function or ``nufeb/volume`` compute.

.. note::
   We are explicitly summing the individual species counts rather than using the ``all`` meta-group. This avoids issues with inadvertently counting particles which should be excluded, like EPS; an error which can easily occur when modifying runs to include such particles after initial prototyping.

.. seealso::
   Details on how LAMMPS handles output, including variables
        `LAMMPS HowTo: output <https://docs.lammps.org/Howto_output.html>`_

   Computing biomass volume in NUFEB
        :doc:`compute_volume`

