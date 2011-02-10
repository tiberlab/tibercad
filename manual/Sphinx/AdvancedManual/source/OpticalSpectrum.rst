.. _OpticalSpectrum:

Simulation opticalspectrum
==========================

By defining the model **opticalspectrum**, optical matrix elements are used to calculate
the associated (emission) spectrum with a k-space integration. In **Solver** section::

  opticalspectrum
    {
     k_space_dimension = 2
     k-space_basis = true
     k1 = (0, 0, 0.1)
     k2 = (0, 0.1, 0)
     refine_fraction = 0.30
     relative_accuracy = 0.01
     refine_k_space = true
     number_of_nodes = (2, 2)
     wedge = quarter
     optical_matr_elem_model = opticskp
     polarization = (0, 0, 1)
     Emin = 3.0
     Emax = 5.0
     dE = 0.001
    }

The parameters:

*k_space_dimension* = **1** for 2D simulations, **2** for 1D simulations. *k-space_basis* is
**true** if the k-space is defined by means of k-vectors; if **false**, vectors are expressed in
real space.

If *refine_k_space* = **true**, that is adaptive k-mesh refinement is enabled, 
all the elements whose error is greater than the value (*1-refine fraction*)* (maximum error) are
going to be refined. In this case, 'Error' is just the integrated quantity. The refinement
will end when the *relative_accuracy* is obtained.

*number_of_nodes* = numb. of elements in k mesh, along each direction

*wedge* = half | quarter, to reduce calculation time, by exploiting symmetry.

*optical_matr_elem_model* = name of the *opticskp* model associated

*polarization* = light polarization (vector)

*Emin, Emax, dE* : energy range and step of spectrum calculation.

Output
------

The output variables for optics calculations are:

* **optical_spectrum** : k-space integrated optical emission spectrum, calculated by
**opticalspectrum** model




.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes