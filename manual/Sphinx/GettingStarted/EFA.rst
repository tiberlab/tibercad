..   <marker>

.. _EFAReferenceGuide:


Quantum  EFA  calculations
=================================


In  tiberCAD,  it  is  possible  to  perform quantum  calculations in  the  framework  of  Envelope Function Approximation (EFA):  eigenstates, eigenfunctions and  quantum  density of  a  given system and    dispersion  of  quantum  states  can  be  obtained  by means  of the module:

* **Module**  ``efaschroedinger``


The  optical properties  are  calculated by  the  module 

* **Module**  ``opticskp``



Module efaschroedinger
-----------------------
 

The  ``efaschroedinger``  simulation tool of tiberCAD  is developed
in order to solve a single-particle    :math:`Schr\ddot{o}dinger`  equation for electrons and holes in a semiconductor crystal. 
This problem is an eigenvalue problem that is treated as a generalized
complex eigenvalue problem

..  math::
    :nowrap:
    :label:

	
    \begin{equation}
	 H\psi  = E S \psi,
    \end{equation}

where H and S are the Hamiltonian and S-matrix, respectively.
The EFA calculations are performed by the **Module** ``efaschroedinger``.


A typical example is the following::


  Module efaschroedinger
  {
    particle = el
    poisson_model_name = dd
    strain_model_name = strain # 
    name = quantum_el
    regions = quantum
    plot = (EigenFunctions, EigenEnergy, EnergyLevels, QuantumDensity)
    
    Solver
    {
      number_of_eigenstates = 10 # 30 
     }
    Physics
    {
      model = single_band 
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
   
    model = kp   #  k.p  for  valence  band
    kp_model = 6x6
   }


Calculation of quantum dispersion
^^^^^^^^^^^^^^^^^^^^^^^

Within  the **Module**  ``efaschroedinger``  it is possible to calculate the dependence of quantum eigenstates on ``k``-vector. Such dependence gives the *quantum state dispersion*. 
To  calculate the  quantum dispersion we  need  to  define the  block *Dispersion*


::

  Dispersion 
  {
    k-path = G-K-M
    number_of_nodes = 10
    k_max = 0.1  
  }


The dispersion is calculated for each of  the quantum states previously  defined in  the  simulation.





Output
^^^^^^^^^^^^^^^^^^^^^^^



The available output variables, to be specified in the plot option, are the following:

  ``EigenEnergy``: 
    Eigen energy in eV

  ``ProbabilityDensity``: 
    square module :math:`|\psi({\bf r})|^2` of the eigenstate wavefunction 

  ``k-space_dispersion`` :
    output of  the quantum dispersion

  ``QuantumDensity``:   
    if  present, quantum  density is  calculated 

  ``Occupation``: 
    probability to find the state occupied. It is calculated assuming Fermi
    distribution and mean electrochemical potential and temperature.







Module opticskp
----------------------


By defining the **Module** ``opticskp`` , calculation of optical properties is enabled; in particular, the optical ``kp`` matrix elements are calculated from the quantum models specified in the **Module**.





The optical spectrum from spontaneous emission is calculated in the following way:

..  math::
    :nowrap:
    :label:


    \begin{equation}
         P(\hbar \omega) = \sum_{i,j} \frac{1}{2\pi^2}  \frac{\omega^2_{ij} e^2 }{m^2 c^3}  
         |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
         \frac{\Gamma/2} {(\hbar \omega_{ij} - \hbar \omega)^2 + (\Gamma/2)^2} d\Omega,
    \end{equation}


where :math:`f_i` and :math:`f_j` are the Fermi distributions and :math:`M_{i,j}` is the optical matrix element between the states  :math:`i` and :math:`j`.



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
By  default, all  the eigenstates calculated in  the electron  and  hole quantum simulations are taken in  account  for the  optical calculations.

By specifying a range of energy values in this way::

  Emin = 3.0
  Emax = 5.0
  dE = 0.001

the emission optical spectrum for **k=0** is calculated.









Integrated spectrum
^^^^^^^^^^^^^^^^^^^^^^^


For  1D and  2D  calculations,  it  is  possible  to  perform  an  integration  of  the  optical  spectrum  in  *k-space*, by specifying *optical_spectrum* in  the  ``plot``  statement ::

  plot = (optical_spectrum)


In  this  case,  one  has  to  define a  **k-integration**  block inside  **Module** ``opticskp``, in this  way ::

  k-integration
  {
   
     k_max = 0.05     
     number_of_elements = (5,5)
     quadrature_type = gaussian 
     quadrature_order = third 

     refine_k_space = false 
     refine_fraction = 0.5
     relative_accuracy = 0.001 

  }


To  run  an  opticskp  sumulation,  one  needs  to  define it in the  solve  statement,  in  this  way::





   Simulation
   { 
      solve = (strain, dd,optics) 
   }


where  *optics*  is  the  name  of  the  defined ``opticskp``  simulation.



Output
^^^^^^^^^^^^^^^^^^^^^^^


The output variables for optics calculations are:

 ``optical_spectrum_k_0`` : 
    optical emission spectrum for *k=0*.

 ``optical_spectrum`` : 
    optical emission spectrum integrated in  k-space.




.. rubric:: Footnotes


..   </marker>
