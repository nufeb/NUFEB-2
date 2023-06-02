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

        *diffdt* value = time step for chemical processes
        *difftol* value = stopping tolerance for chemical processes
        *diffmax* = max iterations for chemical processes
        *pairdt* = time step for physical processes
        *pairtol* = stopping tolerance for physical processes
        *pair max* = max iterations for physical processes


Examples
""""""""

.. code-block::

   group bac type 1
   grid_style nufeb/chemostat 1 glucose 0.001

   fix f_monod bac nufeb/growth/monod glucose 3.5e-5 growth 4e-4 yield 0.61 decay 2e-5 maintain 1e-5
   fix f_monod bac nufeb/growth/monod glucose 3.5e-5 growth 4e-4 yield 0.61