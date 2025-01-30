.. index:: fix nufeb/T6SS/lysis

fix nufeb/T6SS/lysis command
============================

Syntax
""""""

.. parsed-literal::
    
     fix ID group-ID nufeb/T6SS/lysis sub-ID lysis-rate sub-ratio

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* sub-ID = substrate ID for organic substrate released by lysis
* lysis-rate = rate at which lysed organism loses mass 
* sub-ratio = conversion ratio of lost mass to substrate produced
         
Examples
""""""""

.. code-block::

   # lyse at a rate of 8 h^-1 and 20% of mass lost becomes substrate
   fix lysis_intox het_intoxicated nufeb/T6SS/lysis sub 2e-3 0.2
   
Description
"""""""""""
The atoms defined in *group-ID* will undergo lysis, resulting in a mass loss and the conversion of some of that mass into a substrate.
The affected atoms are considered as heterotrophic bacteria, having a spherical shape (:doc:`atom_style coccus <atom_vec_coccus>`).

The mass loss is calculated by the ``lysis-rate`` parameter, which is in units of inverse time. 

.. math::

   m_{t+1} = m_{t}(1-r_{lys} \Delta_t) 


where:

* :math:`m_{t+1}` is the mass of the atom after the time step
* :math:`m_{t}` is the mass of the atom before the time step
* :math:`r_{lys}` is the lysis rate 
* :math:`\Delta_t` is the duration of the time step 

A proportion of the lost mass is released as a substrate, given the ratio ``sub-ratio``, with substrate specified by ``sub-ID``.

The fix is called at each biological step (see :doc:`run_style nufeb <run_style_nufeb>`)
to update atom and grid attributes.

The value of the ``sub-ID`` must be consistent with the name defined in the
:doc:`grid_style chemostat <grid_style_chemostat>` command.


This fix requires that NUFEB be installed/built by including the ``--enable-t6ss`` argument when installing via ``install.sh`` or by including ``make yes-t6ss`` as part of the build process. 

Restrictions
"""""""""""""
This fix is not compatible with the following commands:

* :doc:`atom_style bacillus <atom_vec_bacillus>`


