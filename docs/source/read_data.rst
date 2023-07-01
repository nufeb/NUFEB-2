.. index:: read_data

read_data* command
===================

Syntax
""""""

.. parsed-literal::

   read_data file keyword args ...

* file = name of data file to read in
* zero or more keyword/arg pairs may be appended
* keyword (see `read_data* <https://docs.lammps.org/read_data.html>`_) = *add* or *offset* or *shift* or *extra/atom/types* or *extra/bond/types* or *extra/angle/types* or *extra/dihedral/types* or *extra/improper/types* or *extra/bond/per/atom* or *extra/angle/per/atom* or *extra/dihedral/per/atom* or *extra/improper/per/atom* or *group* or *nocoeff* or *fix*

Examples
""""""""

.. code-block::

   read_data data.nufeb

.. note::

    This page provide information used for NUFEB microbial simulation only.
    More descriptions regarding file format and command keywords are detailed in `read_data* <https://docs.lammps.org/read_data.html>`_.

Read in a data file containing information NUFEB/LAMMPS needs to run a
simulation.  The file can be ASCII text or a gzipped text file
(detected by a .gz suffix).  This is one of 3 ways to specify initial
atom coordinates; see the `read_restart* <https://docs.lammps.org/read_restart.html>`_ and
`create_atom* <https://docs.lammps.org/create_atom.html>`_ commands for alternative methods.
Also see the explanation of the `-restart command-line switch* <https://docs.lammps.org/Run_options.html>`_
which can convert a restart file to a data file.

----------

Format of a data file
"""""""""""""""""""""

A data file has a header and a body.  The header appears first.  The
first line of the header is always skipped; it typically contains a
description of the file.  Then lines are read one at a time.  Lines
can have a trailing comment starting with '#' that is ignored.  If the
line is blank (only white-space after comment is deleted), it is
skipped.  If the line contains a header keyword, the corresponding
value(s) is read from the line.  If it does not contain a header
keyword, the line begins the body of the file.

The body of the file contains zero or more sections.  The first line
of a section has only a keyword.  This line can have a trailing
comment starting with '#' that is either ignored or can be used to
check for a style match, as described below.  The next line is
skipped.  The remaining lines of the section contain values.  The
number of lines depends on the section keyword as described below.
Zero or more blank lines can be used between sections.  Sections can
appear in any order, with a few exceptions as noted below.

----------

Format of the header of a data file
"""""""""""""""""""""""""""""""""""

These are the recognized header keywords for NUFEB simulation.

* *atoms* = # of atoms in system
* *atom types* = # of atom types in system
* *bacilli* = # of bacilli in system
* *xlo xhi* = simulation box boundaries in x dimension
* *ylo yhi* = simulation box boundaries in y dimension
* *zlo zhi* = simulation box boundaries in z dimension

The initial simulation box size is determined by the lo/hi settings.
The system may be periodic or non-periodic; see the
`boundary* <https://docs.lammps.org/boundary.html>`_.
When the simulation box is created
it is also partitioned into a regular 3d grid of rectangular bricks,
one per processor, based on the number of processors being used and
the settings of the `processors* <https://docs.lammps.org/processors.html>`_.

The "bacilli" setting is only used with
:doc:`atom_style bacillus <atom_vec_bacillus>`
and specify how many of the atoms are
bacilli. See the discussion of bacillusflag and the *Bacilli* section below.

Examples:

.. parsed-literal::

   5 atoms
   2 atom types

   0.0   1e-4    xlo xhi
   0.0   0.5e-4  ylo yhi
   0.0   1e-4    zlo zhi

----------

Format of the body of a data file
"""""""""""""""""""""""""""""""""

These are the section keywords for the body of the file (NUFEB-specific only):

----------

*Atoms* section:

* one line per atom
* line syntax: depends on atom style

.. list-table::

   * - bacillus
     - atom-ID atom-type bacillusflag density x y z
   * - coccus
     - atom-ID atom-type diameter density x y z outer_diameter

An *Atoms* section must appear in the data file if natoms > 0 in the header section.
These are the
line formats for each atom_style (:doc:`atom_style coccus <atom_vec_coccus>` or
:doc:`atom_style bacillus <atom_vec_bacillus>`).

.. parsed-literal::

    * atom-ID = integer ID of atom
    * atom-type = type of atom (1-Ntype)
    * bacillusflag = 1 for bacillus particles
    * density = density of particle (kg/m3)
    * diameter = diameter of coccus or width of bacillus (m)
    * outer_diameter = outer_diameter of coccus
    * x,y,z = coordinates of atom (m)

The *atom-ID* is used to identify the atom throughout the simulation and
in dump files.  Normally, it is a unique value from 1 to Natoms for
each atom.  Unique values larger than Natoms can be used, but they
will cause extra memory to be allocated on each processor, if an atom
map array is used, but not if an atom map hash is used; see the
`atom_modify* <https://docs.lammps.org/atom_modify.html>`_
command for details.  If an atom map is
not used (e.g. an atomic system with no bonds), and you don't care if
unique atom IDs appear in dump files, then the atom-IDs can all be set
to 0.

The *diameter* specifies the size of a finite-size spherical coccus if
the command
:doc:`atom_style coccus <atom_vec_coccus>` is used, or the width
of rod-shaped bacillus if :doc:`atom_style coccus <atom_vec_bacillus>`
is used.

The *bacillusflag* determine
whether the particle is a finite-size bacillus.
Additional attributes must be defined for each bacillus
in the corresponding "Bacilli" section.

The *density* is used in conjunction with the
particle volume to set the mass of each particle as mass = density \*
volume.

Examples:

.. parsed-literal::

    Atoms

        1 1 1.0e-6 150 0.5e-5 0.5e-5 1e-6 1.0e-6
        2 1 1.0e-6 150 1.5e-5 0.5e-5 1e-6 1.0e-6
        3 2 1.0e-6 100 2.5e-5 0.5e-5 1e-6 1.2e-6
        4 2 1.0e-6 100 3.5e-5 0.5e-5 1e-6 1.2e-6
        5 2 1.0e-6 100 4.5e-5 0.5e-5 1e-6 1.2e-6


----------

*Bacilli* section:

* one or more lines per body
* first line syntax: atom-ID ixx iyy izz ixy ixz iyz px py pz diameter

.. parsed-literal::

    * atom-ID = integer ID of atom
    * ixx iyy izz ixy ixz iyz = 6 moments of inertia
    * px py pz = coordinate for one of the two poles
    * diameter = width of bacillus (m)

Keywords *ixx iyy izz ixy ixz iyz* sets the 6 moments of inertia for each bacillus.
The values should be consistent with the current orientation of the bacillus around its center of mass.
The values are with respect to the simulation box XYZ axes,
not with respect to the principal axes of the particle itself.
NUFEB performs the latter calculation internally.
The center-of-mass position of the particle is specified by the x,y,z values
in the Atoms section of the data file.

*px py pz* specify the coordinate for one of the two bacillus poles, i.e, the central point of
the hemisphere. The coordinate of the other pole is calculated internally.

For example, a bacillus whose diameter is 1e-6 (width), length 4e-6 (exclude hemispherical caps), and
located in parallel with z-surface along with x-axis, is
specific as follows:

.. parsed-literal::

    Bacilli

      1 0 0 7.2e-24 0 0 0 2.0e-6 0 0 1e-6