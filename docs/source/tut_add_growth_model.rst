Extend NUFEB with new growth model 
===================================

This tutorial provides guidance on developing your own microbial system with new species types and their catabolic models. 
This requires some coding works to implement their unique growth equations into the NUFEB codebase. 
Before doing, make sure none of the growth model provided in NUFEB (fix nufeb/growth/\*) can be used in your system.

Define system
"""""""""""""

Let's suppose your system contains two types of species named 'spec1' and 'spec2' . 
Both of them utilise nutrient 'nut' for their growth, and the growth of spec2 is further inhibited by antibiotics 'anti'. 
The corresponding ODEs are as follow:

.. math::

  \dot m_{spec1} & =  \mu_{max1} \frac{S_{nut}}{S_{nut} + Ks_{nut}} m_{spec1}
  
  \dot m_{spec2} & =  \mu_{max2} \frac{S_{nut}}{S_{nut} + Ks_{nut}} \frac{Ks_{anti}}{S_{anti} + Ks_{anti}} m_{spec2}
   
  \dot S_{nut} & = - (\frac{1}{Y_{spec1}} \dot m_{spec1} + \frac{1}{Y_{spec2}} \dot m_{spec2}) 
  
  \dot S_{anti} & = 0
  
where

* :math:`\mu_{max1}` and :math:`\mu_{max1}` are the maximum growth rates of spec1 and spec2, respectively
* :math:`S_{nut}` and :math:`S_{anti}` are the concentrations of the nutrient and antibiotics, respectively
* :math:`Ks_{nut}` and :math:`Ks_{anti}` are the half-velocity constants of the nutrient and antibiotics, respectively
* :math:`Y_{spec1}` and :math:`Y_{spec2}` are the yield coefficients of the two species

Prepare input files
"""""""""""""""""""

Having the system model defined, the next step is to prepare an inputscript with the new species, nutrients and their parameters.
It is always good to use an existing NUFEB input file as template. Here we use the files from examples/simple-growth 


