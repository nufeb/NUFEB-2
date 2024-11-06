.. index:: compute nufeb/volume

compute nufeb/volume command
=============================

Syntax
------

```compute ID group-ID nufeb/volume```

* ID = User-assigned name for this compute
* group-ID = ID of the group of atoms to apply this compute to

Examples
--------

::

    # Compute the volume of all bacteria 
    compute biomass_vol all nufeb/volume

.. seealso::
    Used in this recipe:
       :doc:`cook_percent_biomass`

Description
-----------

The ``vol`` compute simple returns the volume of all atoms matching the given ``group-ID`` and respects the magic LAMMPS ``all`` group. This compute works for both coccus and bacillus cell morphologies.
