.. index:: grid_modify

grid_modify command
==========================


Syntax
""""""

.. parsed-literal::

    grid_modify set subID xboundary yboundary zboundary S\ :sub:`domain` S\ :sub:`bulk`
    
* subID = ID of the substrate to apply the command to
* xboundary, yboundary, zboundary = *p* or *n* or *d*, two letters 
* S\ :sub:`domain` = initial substrate concentration in the simulation box.
* S\ :sub:`bulk` = initial substrate concentration in bulk liquid

	.. parsed-literal::
	
	    p is periodic
	    n is no-flux 
	    d is Dirichlet
	    
Examples
""""""""

.. code-block:: 

   grid_modify set o2 pp pp nd 1e-4 1e-4
   grid_modify set nh4 pp pp nd 1e-3 1e-3
   grid_modify set glucose nn nn nn 1e-4 0.0
   
Description
""""""""""""""

Set initial values of grid attributes w.r.t, the substrate *subID* defined in grid_style commands
(:doc:`grid_style nufeb/simple <grid_style_simple>` or :doc:`grid_style nufeb/chemostat <grid_style_chemostat>`).

*xboundary*, *yboundary* and *zboundary* set the style of substrate diffusion boundary conditions
in each dimension of the simulation box. 
For each keyword, two letters assigns the first style to the lower face
and the second style to the upper face. 

The style *p* means periodic boundary condition, 
which models a cyclic situation of the substrate concentration flux across the boundary surface. 
The style must be applied to both faces of a dimension. 
For example, *pp* is valid while *pn* or *dp* is invalid.

The style *n* means Neumann or no-flux boundary condition.
It intuitively means that the surface is insulated, and no substrate enters or leaves.

The style *d* means Dirichlet boundary condition,
which intuitively means that the boundary is connecting to a well-mixed bulk liquid, and
the substrate concentration at the boundary is kept fixed.

The combinations of boundary conditions allow to model different simulation scenarios.
For example, 1) *nn nn nn* often model a system that microbes grow in a complete insulated micro well (closed system).
2) *pp pp nd* models a representative volume of the full-scale bioreactor, 
where the (modelled) microbes grows on the surface of the bioreactor. 
3) *dd dd dd* can be applied when modelling a microbial cluster (floc) floating in water.

S\ :sub:`domain` defines the initial substrate concentration inside the simulation box.
S\ :sub:`bulk` defines the initial substrate concentrations in the Dirichlet *d* boundary, i.e, 
initial bulk liquid concentration.



