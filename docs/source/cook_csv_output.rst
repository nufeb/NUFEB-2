Saving simulation output to a CSV file
======================================


Problem
-------

You want to save the simulation output to a CSV file. This can include one or more computations specified in the inputfile.

Solution
--------

Near the end of the inputfile, uses the LAMMPS ``print`` fix to write one or more pre-specified variables to a file.
In the example below, the variables ``step``, ``simvol``, ``fillpct``, ``bug1_relab``, and ``bug2_relab`` are assumed to be valid variables.

.. code-block::

    # Record the cell volume data to a csv file                                     
    fix volcsv all print 100 "${step},${simvol},${fillpct},${bug1_relab},${bug2_relab}" &
        screen no &                                                                 
        file "cell_rel_volumes.csv" &                                               
        title "step,simulation_volume,fill_percent,bug1_relab,bug2_relab"

In the first line, we declare that the fix will be run every 100 timesteps and give the order in which the variables will appear on each data line.

On the second line, we specify to not send the output to the screen.

On the third line, we declare the name of the output file.

The final line provides the header row for the CSV file.

Discussion
----------

During execution of the inputscript, many values can be computed for output at each timestep using computes and variables. Often these values are sent to the screen for immediate output via the LAMMPS ``thermo`` command, but sometimes it would be preferable to save them into a well-formatted file for later analysis. Taking advantage of the LAMMPS ``print`` fix is preferable than trying to redirect the output to a file.

.. note::
   The example uses the ``&`` line continuation symbol to aid in readability.

.. seealso::
   Details on how LAMMPS handles output
        `LAMMPS HowTo: output <https://docs.lammps.org/Howto_output.html>`_

