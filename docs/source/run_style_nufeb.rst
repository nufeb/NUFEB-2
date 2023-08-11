.. index:: run_style nufeb

run_style nufeb command
========================

Syntax
""""""

.. parsed-literal::

    run_style nufeb keyword values ...

* one or more keyword/value pairs may be listed
* group-ID = ID of the group atoms to apply the fix to

    .. parsed-literal::

        *diffdt* value = time step for diffusion and chemical processes (default: 1.0e-3 s)
        *difftol* value = stopping tolerance for diffusion and chemical processes (1.0e-6 kg/m3)
        *diffmax* = maximum # of iterations for diffusion and chemical processes (default: -1)
        *pairdt* = time step for physical processes (default: 1.0e-8 s)
        *pairtol* = stopping pressure tolerance for physical processes (default: 1.0 N/m2)
        *pairmax* = maximum # of iterations for physical processes (default: -1)
        *profile* = *file_name*, print performance info to file
        *screen* = *yes* or *no*, print additional diffusion and pressure information to screen (default: yes)
        *initdiff* =  *yes* or *no*, solve diffusion during initialisation (default: yes)


Examples
""""""""

.. code-block::

   #--- examples/biofilm-heterotroph ---#

   run_style nufeb diffdt 1e-4 difftol 1e-6 pairdt 1e-2 pairtol 1 pairmax 1000 diffmax 5000


Description
""""""""""""""

Choose the style of time integrator used for a NUFEB simulation.
The procedure of the integration is shown below:

.. parsed-literal::

    Run for N `biological steps <https://docs.lammps.org/run.html>`_ and for each iteration, do:
        1. change to `biological timestep <https://docs.lammps.org/timestep.html>`_
        2. solve :doc:`biological processes <list_biology>`
        3. change to physical timestep *pairdt*
        4. solve :doc:`physical processes <list_physics>` until reaching *pairtol* or *pairmax*
        5. change to `biological timestep <https://docs.lammps.org/timestep.html>`_
        6. solve :doc:`post-physical processes <list_post_physics>`
        7. change to diffusion timestep *diffdt*
        8. solve :doc:`diffusion and chemical processes <list_chemistry>` until reaching *difftol* or *diffmax*
        9. change to `biological timestep <https://docs.lammps.org/timestep.html>`_
        10. solve :doc:`reactor processes <list_reactor>`
        11. :doc:`output <list_post>` results


The model processes are operated sequentially and they are on different timescales.
The coupling between the multiple timescales relies on the pseudo steady-state approximation and the frozen state.
For example, when a threshold pressure (step 4) or a steady state substrate concentration  (step 8) is reached,
the system is assumed to remain unchanged (frozen state) until the next biological step.

The biological processes (step 2) update
the biological attributes of the particles (mass, diameter, type, etc).
The corresponding timestep and the maximum number of biological steps are defined by the LAMMPS commands
`timestep* <https://docs.lammps.org/timestep.html>`_ and
`run* <https://docs.lammps.org/run.html>`_, respectively.
Thus, the total simulation time is calculated as *bio_timestep \* bio_steps*.

The physical processes (step 4) apply the Verlet algorithm to update the physical attributes of
the particles (velocity, force, position, etc).
The algorithm involves a sub-loop using the timestep of *pairdt*.
The keywords *pairtol* and *pairmax* determine the criterion to stop the Verlet integration.
*pairtol* defines a threshold pressure, assuming that the system reaches a mechanical relaxed state
when the total pressure is below the value.
*pairmax* defines the maximum number of physical iterations.
The procedure moves to step 5 when the Verlet sub-loop reaches the specified maximum number,
regardless of the current system pressure.

The post-physical processes (step 6) update
particle or grid attributes that need to be considered after solving the physical processes.
A typical example is :doc:`fix nufeb/diffusion_coeff <fix_diffusion_coeff>`.
The command updates diffusion coefficient at each grid based on its occupied particles.
Therefore, a physically relaxed system is required in order to accurately calculate the values.

The chemical processes (step 8) are solved on a Marker-And-Cell (MAC)
uniform grid (see :doc:`grid style <grid_style_chemostat>`)
using a :doc:`diffusion solver <fix_diffusion>`.
The chemical attributes of each grid (concentration, pH, reaction rate, etc)
are updated during each diffusion iteration.
The keywords *difftol* and *diffmax* determine the criterion for stopping the diffusion
and all other chemical processes.
*pairtol* defines a threshold substrate concentration
used to evaluate whether the system reaches a steady state substrate concentration.
*diffmax* is the maximum number of diffusion iteration.

The reactor processes (step 10) update the attributes of the large-scale system (reactor) to which the
microscale IbM simulation box is connected (e.g, via dirichlet boundary condition).
Example attributes include bulk solute concentration (:doc:`nufeb/reactor/solute_balance <fix_reactor_solute_balance>`),
bulk gas concentration (:doc:`nufeb/reactor/gas_balance <fix_reactor_gas_balance>`),
and boundary layer position (:doc:`nufeb/boundary_layer <fix_boundary_layer>`).

The *profile* keyword outputs the performance of each IbM process module to a file.
When the *screen* keyword is enabled, additional diffusion and pressure information is displayed in the terminal after
each biological step.
When the *initdiff* keyword is activated, the diffusion solver will be triggered during the simulation initialisation stage.
This allows for the updating of substrate concentration before addressing the biological processes.