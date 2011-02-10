.. _Opticskp:

Simulation opticskp
===================

By defining the **opticskp** model, calculation of optical properties is enabled; in particular, 
the optical kp matrix elements are calculated from the quantum models specified in

**Solver** section::

  opticskp
    {
     initial_state_model = QW1_electrons # quantum_el
     final_state_model = QW1_holes # quantum_hl
     initial_eigenstates = (0, 19)
     final_eigenstates = (0, 19)
    }
  
Here, *initial_state_model* and *final_state_model* are, respectively, the quantum models
(**efaschroedinger** model) associated respectively to the initial state of optical transition
(e.g. electron), and to the final state of optical transition (e.g. hole). *initial_eigenstates*
and *final_eigenstates* refer to the range of eigenstates to be taken in account for optical
calculations. By specifying, in **Solver section**, a range of energy values in this way::

  Emin = 3.0
  Emax = 5.0
  dE = 0.001
  
the emission optical spectrum for **k=0** is calculated. The spectrum is calculated in the
following way:

.. math:.
   :nowrap:
   :label:
   
   missing Source equation   (8.1)
   
where :math:`f_i` and :math:`f_j` are the Fermi distributions.


Output
------

The output variables for optics calculations are:

* **optical_spectrum_k_0** : optical emission spectrum for **k=0** , calculated through **opticskp** model





.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes