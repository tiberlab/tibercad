.. _QuantumDensity:

Quantum density
===============

This is an example of Quantum Density::

  dens_el
    {
     k_space_dimension = 2
     k-space_basis = true
     k1 = (0, 0, 0.1)
     k2 = (0, 0.1, 0)
     number_of_nodes = (4,4)
     wedge = quarter
     refine_fraction = 0.20
     relative_accuracy = 0.01
     refine_k_space = true
     uniform_refinement = false
     mesh_order = FIRST
     quantum_simulation = quantum_el
     degeneracy = 2
     initial_eigenstates_number = 10
     analitic = false
    }
  
* **quantum simulation** name of the Schroedinger simulation

.. index:: single:Degeneracy

* **:term:`degeneracy`** degeneracy of the quantum state

* **initial eigenstates number** initial number of eigenstates for the Schroedinger
equation

* **analytic** = { true | false } If true then the density is calculated analytically or
numerically. 

Analitical calculation of density is done in the following way. For each eigenstate
we calculate the effective mass assuming quadratic dispersion. Then the charge
density is calculate as follows:

.. math::
   :nowrap:
   :label:
   
   missing equation source (11.2)
w
here :math:`\rho_{1D}` and :math:`\rho_{2D}` are the 1D and 2D charge densities; *m* is the averaged mass
(the mass is different for each quantized state and is position independent); *g* is
the degeneracy of the states. The :math:`+` sign is for electrons, the :math:`-`  sign is for holes.
Numerical calculation is done by the following formula::

.. math::
   :nowrap:
   :label:
   
   missing equation source(11.3)
   
The integration is performed on a mesh in the k-space.
Output
------

..  index:: single:quantum density

The output parameter is **quantum density** .


.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes