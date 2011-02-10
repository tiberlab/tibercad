.. _QuantumDispersion:

Quantum dispersion
==================

There is a possibility to calculate the dependence of quantum eigenstates on k-vector.
Such dependence is called *dispersion*. The simulation name is **quantumdispersion**.

Example::

  model quantumdispersion
    {
     options
       {
        simulation_name = dispersion1D_el
        physical_regions = all
       }
    }
  
Solver options
--------------

The dispersion of quantum states is calculate at k-points that are nodes of the mesh in k-space:

* **quantum_simulation** name of the Schroedinger equation simulation

* **min_eigenvalue_number**, **max eigenvalue_number** : the dispersion is calculated for the states number i , where
 
max_eigenvalue_number :math:`\geq  i \geq` min_eigenvalue_number
 
The rest of the parameters (wedge, k space dimension, etc...) define the k-space::

  dispersion1D_el
    {
     quantum_simulation = quantum_el
     min_eigenvalue_number = 0
     max_eigenvalue_number = 5
     wedge = half
     k_space_dimension = 1
     k1 = (0, 0.1, 0)
     number_of_nodes = (10)
     output_format = grace
    }

Output
------

The output variable name is **k-space_dispersion**. The output format for the dispersion
can be controlled independently of the general specification in the Simulation section
by redefining the *output_format* keyword.


.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes