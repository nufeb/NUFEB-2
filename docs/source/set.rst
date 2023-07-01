.. index:: set

set* command
=============

Syntax
""""""

.. parsed-literal::

   set style ID keyword values ...

* style = *atom* or *type* or *group* or *region*
* ID = depends on style

.. parsed-literal::

       for style = *atom*, ID = a range of atom IDs
       for style = *type*, ID = a range of numeric types or a single type label
       for style = *group*, ID = a group ID
       for style = *region*, ID = a region ID

* one or more keyword/value pairs may be appended
* keyword = *x* or *y* or *z* or *diameter* or *density* or *mass* or *outer_mass* or *outer_diameter* or *outer_density* or *biomass* or *bacillus/inertia*
  or *bacillus/pole/random* or *bacillus/length*
  *type* or *type/fraction* or *type/ratio* or *type/subset*

.. parsed-literal::

    *x*,\ *y*,\ *z* value = atom coordinate (m)
    *diameter* value = diameter of spherical coccus, or width of rod bacillus (m)
    *density* value = atom density (kg/m3)
    *mass* value = wet mass (kg)
    *outer_mass* value =  EPS shell mass (exclude inner mass)
    *outer_diameter* value  = outer diameter (diameter + EPS shell depth)
    *outer_density* value = EPS shell density
    *biomass* value = dry mass ratio
    *bacillus/length* value = length of bacillus
    *bacillus/inertia* value = ixx iyy izz ixy ixz iyz
     ixx iyy izz ixy ixz iyz = 6 moments of inertia
    *bacillus/pole/random* value = angle seed
     angle = orientation of line segment with respect to one of the 7 directions: *x* or *y* or *z* or *xy* or *yz* or *xz* or *xyz*
     seed = random # seed (positive integer) for line segment orientations
    *type* value = numeric atom type or type label
     value can be an atom-style variable

Examples
""""""""

.. code-block::

    set type 1 z 1e-6
    set type 2 diameter 1e-6
    set group HET outer_diameter 1.2e-6
    set group AOB density 290
    set group AOB NOB biomass 0.6
    set group Ecoli bacillus/length 3e-6
    set group Ecoli diameter 0.8e-6
    set group Ecoli bacillus/inertia 0 0 9e-23 0 0 0
    set group Ecoli bacillus/pole/random xy 1234

.. note::

    This page describes properties used for :doc:`atom_style coccus <atom_vec_coccus>`
    and :doc:`atom_style bacillus <atom_vec_bacillus>` only.
    Properties for other LAMMPS atom types are detailed in `set* <https://docs.lammps.org/set.html>`_.

Set one or more properties of coccus and bacillus atoms. The command can be useful for overriding the default or initial
values assigned by the  `create_atom* <https://docs.lammps.org/crete_atom.html>`_, :doc:`read_data <read_data>`
or `read_restart* <https://docs.lammps.org/read_restart.html>`_ command.

The style *atom* selects all the atoms in a range of atom IDs.

The style *type* selects all the atoms in a range of types or type
labels.  The style *type* selects atoms in one of two ways.  A range
of numeric atom types can be specified.  Or a single atom type label
can be specified, e.g. "C".

In each of the range cases, the range can be specified as a single
numeric value, or a wildcard asterisk can be used to specify a range
of values.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  For
example, for the style *type*, if N = the number of atom types, then
an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).  For all the styles except
*mol*, the lowest value for the wildcard is 1; for *mol* it is 0.

The style *group* selects all the atoms in the specified group.  The
style *region* selects all the atoms in the specified geometric
region.  See the  `region* <https://docs.lammps.org/region.html>`_. and  `group* <https://docs.lammps.org/group.html>`_.
for details of how to specify a group or region.

----------

This section describes the keyword options for which properties to
change, for the selected atoms.

Keyword *type* sets the atom type for all selected atoms. A specified
value can be either a numeric atom type or an atom type label. When
using a numeric type, the specified value must be from 1 to ntypes,
where ntypes was set by the `create_box* <https://docs.lammps.org/create_box.html>`_ command or
the *atom types* field in the header of the data file read by the
`read_data* <https://docs.lammps.org/read_data.html>`_ command.

Keywords *x*, *y*, *z* set the coordinates of all selected atoms.

Keyword *diameter* sets the size of the selected atoms. The value must be positive.
Note that this command does not adjust the atom mass.

Keyword *density* also sets the mass of all selected atoms if they have a positive *diameter* value.
Therefore, the *diameter* of atoms must already be defined in order to set the mass using this command.

Keyword *mass* sets the mass of all selected atoms. The value must be positive.

Keyword *outer_diameter* sets the outer diameter of the selected atoms defined as :doc:`atom_style coccus <atom_vec_coccus>` with EPS shell.
The value is the sum of EPS shell depth and the inner diameter, and must be greater or equal than the inner *diameter*.

Keyword *outer-density* sets the EPS shell density of the selected atoms.
It also set the outer mass of EPS shell for atoms if they have a positive outer diameter attribute.

Keyword *outer_mass* sets the mass of the outer EPS shell of the selected atoms.
The value must be greater than or equal to 0.

Keyword *bacillus/length* sets length of the selected atoms defined a :doc:`atom_style bacillus <atom_vec_bacillus>`.
Since bacillus is represented as a cylinder with hemispherical caps.
The length is the distance between the two hemispherical caps (i.e,
the height of the cylinder).

Keyword *bacillus/inertia* sets the 6 moments of inertia for the atoms defined as :doc:`atom_style bacillus <atom_vec_bacillus>`.
The values should be consistent with the current orientation of the bacillus around its center of mass.
The values are with respect to the simulation box XYZ axes,
not with respect to the principal axes of the particle itself.
NUFEB performs the latter calculation internally.
The center-of-mass position of the particle is specified by the x,y,z values above.

Keyword *bacillus/pole/random* sets initial orientation of atoms defined as :doc:`atom_style bacillus <atom_vec_bacillus>`.
The orientation corresponds to the line segment with respect to one of the 7 directions:

* *x* *y* *z* refer to the x-, y-, and z-axes of the box;
* *xy* *yz* *xz* indicate that the line segment is in parallel to the corresponding surface, while the orientation of the third direction is determined randomly based on *seed* value;
* *xyz* indicates a random orientation with respect to the three axes.