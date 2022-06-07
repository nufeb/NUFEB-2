.. index:: fix nufeb/growth/ecoli

fix nufeb/growth/ecoli command
====================================

Syntax
""""""

.. parsed-literal::
    
    fix ID group-ID nufeb/growth/ecoli suc-ID suc-Ks o2-ID o2-Ks co2-ID keyword value ...

* ID = user-assigned name for the fix
* group-ID = ID of the group atoms to apply the fix to
* suc-ID = ID of the substrate (sucrose) for atom growth 
* suc-Ks = half-velocity constant (Ks) for sucrose
* o2-ID = ID of the substrate (oxygen)
* o2-Ks = half-velocity constant (Ks) for oxygen
* co2-ID = ID of the substrate (carbon dioxide)
* zero or more keyword/value pairs may be appended
* keyword = *growth* or *yield* or *decay* or *maintain* 

	.. parsed-literal::
	
	    *growth* value = maximum growth rate 
	    *yield* value = yield coefficient
	    *decay* value = decay rate
	    *maintain* value = maintenance coefficient

Examples
""""""""

.. code-block:: 

   group ecoli type 1
   grid_style nufeb/chemostat 3 suc o2 co2 4e-6
   
   fix f_gecoli aob nufeb/growth/ecoli suc 3.6 o2 0.001 co2 growth 0.00027 yield 0.43 maintain 9.5e-07 decay 2e-05
   
   
Description
"""""""""""
Perform microbial growth to the atoms defined in *group-ID*. 
The affected atoms are considered as *Escherichia coli* (*E.coli*) ..... TODO