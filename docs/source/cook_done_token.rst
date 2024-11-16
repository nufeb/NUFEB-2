Creating a file to indicate the run completed successfully
==========================================================


Problem
-------

You want to create a persistent indicator that a run was successfully completed. One situation where this is useful is when using a workflow manager, such as Snakemake, which uses the filesystem to infer which tasks to perform.

Solution
--------

At the end of the inputfile, use the LAMMPS built-in ``shell`` command, along with ``touch`` to create the file.

::

    # If we get here, we are successful, so create a file to indicate that is so 
    shell touch done.tkn

Discussion
----------

While NUFEB follows the unix style conventions of returning an exit code of 0 when the run completes successfully, you may also want to create some persistent indicator that the run was performed successfully. This will also ease using file-based workflows, such as make and Snakemake. To do this, we can simply tell NUFEB to create an empty file as its last instruction; in effect signalling everything else ran successfully.  The name of the file can be any valid filename.

The LAMMPS ``shell`` command allows you to pass along arguments to the shell, so we can simply use ``touch`` to create the file.

.. seealso::
    Used in as part of the following tutorial:
        :doc:`tut_multi_snakemake`

    Calling shell commands from LAMMPS
        `LAMMPS shell command  <https://docs.lammps.org/shell.html>`_
