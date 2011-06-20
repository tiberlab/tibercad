..   <marker>

.. _EFAReferenceGuide:


Quantum  EFA  calculations
=================================

In  tiberCAD,  it  is  possible  to  perform quantum  calculations in  the  framework  of  Envelope Function Approximation (EFA):  eigenstates, eigenfunctions and  quantum  density of  a  given system and    dispersion  of  quantum  states  can  be  obtained   respectively by means  of the   following Modules:

* Module efaschroedinger

* Module quantumdispersion


The  optical properties  are  calculated by  the  following modules 

* Module opticskp

* Module opticalspectrum


Module efaschroedinger
-----------------------


The EFA calculation  of eigenstates and  eigenfunctions are performed by the **Module** efaschroedinger.  

A typical example is the following::


  Module efaschroedinger
  {
    particle = el
    poisson_model_name = dd
    strain_model_name = strain # 
    name = quantum_el
    regions = quantum
    plot = (ProbabilityDensity, EigenEnergy, QuantumDensity)
    Solver
    {
      number_of_eigenstates = 10 # 30 
     }
    Physics
    {
      model = conduction_band
     }
   }

In  this  example, Schroedinger  equation is  solved for  electrons with a  single band model.
We  calculate  10 eigenstates  by  specifying in  Solver  section::

  number_of_eigenstates = 10  


The  electron quantum density is  calculated too if  the  keyword QuantumDensity is  present in the plot list::

  plot = (ProbabilityDensity, EigenEnergy, QuantumDensity)  


A similar  definition is  needed  for  the  calculation of hole  quantum  states, by  defining:: 

  particle = hl


For  holes, one  can  choose  a  6 bands ``kp``  model in  this  way::

  Physics
  {
    particle = hl
    model = kp   #  k.p  for  valence  band
    kp_model = 6x6
   }





Output
^^^^^^^^^^^^^^^^^^^^^^^



The available output variables, to be specified in the plot option, are the following:

  ``EigenEnergy``: 
    Eigen energy in eV

  ``ProbabilityDensity``: 
    square module :math:`|\psi({\bf r})|^2` of the eigenstate wavefunction 
  ``QuantumDensity``:   
    if  present, quantum  density is  calculated 

  ``Occupation``: 
    probability to find the state occupied. It is calculated assuming Fermi
    distribution and mean electrochemical potential and temperature:


Module quantumdispersion
-----------------------





With the Module quantumdispersion it is possible to calculate the dependence of quantum eigenstates 
on ``k``-vector. Such dependence gives the *quantum state dispersion* . The
simulation name is ``quantumdispersion`` .

::

  Module quantumdispersion
  {
    simulation_name = dispersion1D_el
    regions = all
    quantum_simulation = quantum_el
    min_eigenvalue_number = 0
    max_eigenvalue_number = 5
    wedge = half
    k_space_dimension = 1
    k1 = (0, 0.1, 0)
    number_of_nodes = (10)
   output_format = grace
    plot = k-space_dispersion
  }

The dispersion of quantum states is calculated at k-points that are nodes of the mesh
in k-space.

The main parameters are:

*  ``quantum simulation`` : 
    name of the efaschroedinger simulation.

*  ``min eigenvalue number`` , ``max eigenvalue number`` : 
    the dispersion is calculated for the states number *i* , where 
    
             ``max_eigenvalue_number`` :math:`\ge i \ge`  ``min_eigenvalue_number``

The rest of the parameters (wedge, k space dimension, etc...) define the k-space.



Output
^^^^^^^^^^^^^^^^^^^^^^^


The output variable name is ``k-space_dispersion`` . The output format for the dispersion
can be controlled independently of the general specification in the ``Simulation`` section
by redefining the ``output_format`` keyword.



Module opticskp
----------------------


By defining the Module ``opticskp`` , calculation of optical properties is enabled; in particular, 
the optical ``kp`` matrix elements are calculated from the quantum models specified in the Module.

The optical spectrum from spontaneous emission is calculated in the following way:

where :math:`f_i` and :math:`f_j` are the Fermi distributions.

::

  Module opticskp
  {
    name = optics
    regions = quantum
    plot = (optical_spectrum_k_0 )
    initial_state_model = quantum_el
    final_state_model = quantum_hl
    initial_eigenstates = (0, 9)
    final_eigenstates = (0, 15)
    polarization = (0, 0, 1)
    Emin = 2.8
    Emax = 3.6
    dE = 0.001
  }

Here, ``initial_state_model`` and ``final_state_model`` are, respectively, the quantum simulations 
( ``efaschroedinger`` module) associated respectively to the initial state of optical
transition (e.g. electron), and to the final state of optical transition (e.g. hole). 

``initial_eigenstates`` and ``final_eigenstates`` refer to the range of eigenstates to be taken in
account for optical calculations.

By specifying a range of energy values in this way::

  Emin = 3.0
  Emax = 5.0
  dE = 0.001

the emission optical spectrum for **k=0** is calculated.

Output
^^^^^^^^^^^^^^^^^^^^^^^


The output variables for optics calculations are:

*  ``optical_spectrum_k_0`` : optical emission spectrum for *k=0*.



Module opticalspectrum
---------------------------



By defining the Module ``opticalspectrum`` , optical matrix elements are used to calculate
the associated (emission) spectrum with a k-space integration.

::

  Module opticalspectrum
  {
    k_space_dimension = 2
    k-space_basis = true
    k1 = (0, 0, 0.1)
    k2 = (0, 0.1, 0)
     
    number_of_nodes = (2, 2)
    wedge = quarter
    plot = (optical_spectrum)
    optical_matr_elem_model = opticskp
    polarization = (0, 0, 1)
    Emin = 3.0
    Emax = 5.0
    dE = 0.001
  }

The parameters are the following

*  ``k_space_dimension`` :
   Options are **1** for 2D simulations, **2** for 1D simulations. 

*  ``k-space basis`` : 
   if **true** then the k-space is defined by means of k-vectors; 
   if **false** , vectors are expressed in real space.

*  ``number_of_nodes`` :
   numb. of elements in k mesh, along each direction

*  ``wedge`` :
   half | quarter, to reduce calculation time, by exploiting symmetry.

*  ``optical_matr_elem_model`` :
   name of the *opticskp* model associated

*  ``polarization`` :
   light polarization (vector)

*  ``Emin, Emax, dE`` : 
   energy range and step of spectrum calculation.

Output
^^^^^^^^^^^^^^^^^^^^^^^

The output variables for optics calculations are:

*  ``optical_spectrum`` : k-space integrated optical emission spectrum. 


.. rubric:: Footnotes


..   </marker>
