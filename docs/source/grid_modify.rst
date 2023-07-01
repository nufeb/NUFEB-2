.. index:: grid_modify

grid_modify command
==========================


Syntax
""""""

.. parsed-literal::

    grid_modify set subID xboundary yboundary zboundary S-init keyword value
    
* subID = ID of the substrate to apply the command to
* xboundary, yboundary, zboundary = *p* (periodic) or *n* (no-flux) or *d* (Dirichlet), two letters
* S-init = initial substrate concentration (kg/m3)
* zero or more keyword/value pairs may be appended
* keyword = *mw* or *bulk*

	.. parsed-literal::

	    *bulk* value = initial substrate concentration in bulk liquid (kg/m3)
	    *mw* value = molecular weight (g/mol)
	    
Examples
""""""""

.. code-block:: 

   grid_modify set o2  pp pp nd 1e-4
   grid_modify set nh4 pp pp nd 0.0 bulk 1e-3
   grid_modify set nh3 pp pp nd 1.7e-3 mw 17.031
   
Description
""""""""""""""

Set initial values of grid and substrate attributes defined in grid_style commands
(:doc:`grid_style nufeb/chemostat <grid_style_chemostat>`).

The keywords *xboundary*, *yboundary* and *zboundary* determine the style of diffusion boundary conditions
in each dimension of the simulation box. 
Each keyword consists of two letters assign the first style to the lower surface
and the second style to the upper surface.

The style *p* represents periodic boundary conditions,
which models a cyclic situation of that the substrate flux across the boundary surface.
The style must be applied to both surfaces of a dimension.
For example, *pp* is valid combination, while *pn* or *dp* are invalid.

The style *n* represents Neumann or no-flux boundary conditions.
It means that the surface is insulated, and no substrate enters or leaves.

The style *d* represents Dirichlet boundary conditions,
which assumed that the boundary is connected to a large scale reactor where the
substrate concentration is fixed and well-mixed.

The combinations of boundary conditions allow to model various simulation scenarios.

For example,

*  *nn nn nn* often used to model a system where microbes grow in a complete insulated micro well (closed system).
*  *pp pp nd* models a representative volume of the full-scale bioreactor, where microbes grows on the bottom surface of the bioreactor.
*  *dd dd dd* can be applied to model a microbial cluster (floc) floating in water.

The keyword *S-init* assigns the initial substrate concentration at each grid cell in the simulation box.

The keyword *bulk* determines the initial substrate concentrations in the Dirichlet *d* boundary, i.e,
initial bulk concentration. If the keyword is not given, the concentration is same as *S-init*

The *mw* keyword defines molecular weight of the substrate. Some commands,
such as :doc:`fix nufeb/growth/energy <fix_growth_energy>`, uses the unit of moles (mol), which requires the molecular weight
for the unit conversion.

