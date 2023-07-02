.. index:: fix nufeb/growth/eps

fix nufeb/growth/eps command
============================

Syntax
""""""

.. parsed-literal::

     fix ID group-ID nufeb/growth/eps sub-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = ID for organic substrate
* zero or more keyword/value pairs may be appended
* keyword = *growth* or *yield* or *decay*

	.. parsed-literal::

	    *decay* value = decay rate

Examples
""""""""

.. code-block::

   group EPS type 1
   grid_style nufeb/chemostat 1 sub

   fix f_geps EPS nufeb/growth/eps sub decay 2e-6

Description
""""""""""""""

Perform EPS (extracellular polymeric substances) decay to the atoms defined in *group-ID*.
The affected atoms are considered as EPS, which have a
spherical shape (see :doc:`atom_style coccus <atom_vec_coccus>`).
The model assumes that EPS dissolves into organic substrate in a constant *decay* rate.

The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom attributes.
The value of the organic substrate ID keyword *sub-ID* must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.
The following forward Euler method is implemented to update the mass
(*m*) of each atom in the group:

.. math::
  m' & = m - b_{decay} \cdot \Delta t

The new mass is then used to update the diameter.

Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`