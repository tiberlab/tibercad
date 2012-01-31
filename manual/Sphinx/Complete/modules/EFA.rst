..   <marker>

.. _EFATheory:


Quantum  EFA  calculations
=================================


In  tiberCAD,  it  is  possible  to  perform quantum  calculations in  the  framework  of  Envelope Function Approximation (EFA):  eigenstates, eigenfunctions and  quantum  density of  a  given system and    dispersion  of  quantum  states  can  be  obtained    by means  of the module:

* Module efaschroedinger


The  optical properties  are  calculated by  the  module 

* Module opticskp








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
A typical example is the following  ::

  Module efaschroedinger
    {
     particle = el
     poisson_model_name = dd
     strain_model_name = strain # macrostrain
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


In  this  example, :math:`Schr\ddot{o}dinger`  equation is  solved for  electrons ::

  particle = el 


with a  single band model ::

 model = single_band


We  calculate  10 eigenstates  by  specifying::

  number_of_eigenstates = 10  




The  electron quantum density is  calculated by default if  the  keyword ``QuantumDensity`` is  present in the plot list::

  plot = (ProbabilityDensity, EigenEnergy, QuantumDensity)  


A similar  definition may be  used  in a **Module** ``efaschroedinger`` for  the  calculation of hole  quantum  states, by  defining :: 

  particle = hl


For  holes, one  can  choose, for example,  a  6 bands :math:`{\bf k \cdot p}`  model in  this  way::

  Physics
  {
    model = kp   #  k.p  for  valence  band
    kp_model = 6x6
  }



Module options
^^^^^^^^^^^^^^^^^^^^^^^



The following options influence the behaviour of the Module efaschroedinger:

 ``particle`` : string  
    defines for which particle (electron or hole) :math:`Schr\ddot{o}dinger` equation is solved. 
    Possible values are el and hl. A different Module efaschroedinger has to be
    defined for each particle to be solved.

 ``poisson_model_name`` : string
    defines the name of the simulation (e.g.  driftdiffusion) that can provide electric potential

 ``strain_model_name`` : string 
    defines the name of the simulation (e.g. elasticity) that can provide elastic strain

 ``regions`` : string 
    defines the regions associated to this EFA simulation


Solver section
^^^^^^^^^^^^^^^^^^^^^^^


The Solver section of the Module efaschroedinger contains the following options:

 ``number_of_eigenstates`` : integer 
    defines the number of eigenvalues and eigenfunctions to be found.

 ``Dirichlet_bc_everywhere`` : boolean 
    if *true* (default value), Dirichlet boundary conditions are imposed over all the boundaries of the simulation region

 ``solver`` : string 
    defines the solver for the eigenvalue problem, possible values are: *arnoldi*, *lapack*, *krylovshur*. The default value is **krylovshur** . 

In the case of the *lapack* solver all the eigenvalues are computed. In the case of *arnoldi* or *krylovshur*
solver it is necessary to specify which and how many eigenvalues have to be computed. 
The idea is that the iterative solver calculates several eigenvalues that are
close to a specific number, referred to as the *guess*.


 ``max_iteration_number`` : integer 
    maximum number of iteration, used as a stop condition

 ``eigen_solver_tolerance`` : double 
    numerical eigensolver tolerance used as a convergence criteria

 ``guess`` : double 
     the algorithm try  to  find  the closest eigenvalues to this absolute  value  of  energy (eV). If not defined, then by  default it is  calculated internally from the band edges.

 ``ksp_type`` : string 
    Krylov subspace method type: *bcgsl*, *gmres*, *cg*

 ``pc_type`` : string 
    preconditioner type: *cholesky*, *jacobi*, *ilu* , *composite*.


Physics section
^^^^^^^^^^^^^^^^^^^^^^^


 ``model`` : string 
    possible values are *single_band* , for single  band model ( :math:`\Gamma` point) ;  *kp* for :math:`{\bf k \cdot p}` model

 ``kp_model`` : string 
    possible values are *6x6*, *8x8*. 


In case *single_band* model is  applied to electrons, the  relevant  mass  is read from  the  material file. In case  of  holes, in  *Physics*  section the  following options have  to be  defined:

 ``mass`` : double(units of :math:`m_0`)
    effective  mass to  be  used  in *single_band* calculations for  holes 

 ``band_edge`` : double
    optional override  of  band edge  energy  defined in the material file (:math:`E_v`)



Quantum density calculation
^^^^^^^^^^^^^^^^^^^^^^^



The particle (electron,hole) **quantum density** is  calculated  by default if  the  keyword ``QuantumDensity`` is  present in the plot list::

  plot = (ProbabilityDensity, EigenEnergy, QuantumDensity)  


 
The calculation  may be obtained with an analytical or a numerical calculation.

By default, an analytical calculation of density is performed in the following way. For each eigenstate we calculate the effective mass assuming quadratic dispersion. Then the charge density is calculated  as follows:

..  math::
    :nowrap:
    :label:


    \begin{align}
    \rho_{1D}({\bf r}) &= g\frac{mkT}{2 \pi \hbar^2}  |\psi({\bf r})|^2 \ln \left(1 + 
    \exp \left( \pm \frac{\mu - E}{kT} \right) \right) \\      
    \rho_{2D}({\bf r}) &= g|\psi({\bf r})|^2 \frac{1}{2} \sqrt{\left( \frac{mkT}{2\pi\hbar^2} \right)}  F_{-1/2}\left( \pm \frac{\mu - E}{kT}         \right),       
    \end{align}


where :math:`\rho_{1D}` and :math:`\rho_{2D}` are the 1D and 2D charge densities; :math:`m` is the averaged mass (the mass
is different for each quantized state and is position independent); :math:`g` is the degeneracy of
the states. The :math:`+` sign is for electrons, the :math:`-` sign is for holes.



Alternatively, numerical calculation may be performed according  the following formula:

..  math::
    :nowrap:
    :label:

    \begin{equation}
    \rho({\bf r}) = \sum_n \frac {1}{(2\pi)^d} \int_{BZ} |\psi_{\bf k_{\|}}|({\bf r})|^2 \frac{1}{1+\exp 
    \left(\pm \frac{E-\mu}{kT} \right)} d{\bf k_{\|}}
    \end{equation}


The integration is performed on a mesh in the k-space.
To perform a numerical calculation,  a  block ``quantumdensity`` has  to  be  defined inside the **Module** ``efaschroedinger``.

:: 

  quantumdensity
    {
     analytic = false
     k_max = 0.05     
     number_of_elements = (5,5)
    }

The available options are:



  ``analytic`` : bool  
      If true then the density is calculated analytically, otherwise numerically.
  ``number_of_elements`` :
      number of elements in the mesh along each  direction





Calculation of quantum dispersion
^^^^^^^^^^^^^^^^^^^^^^^


Within  the **Module**  ``efaschroedinger``  it is possible to calculate the dependence of quantum eigenstates 
on ``k``-vector. Such dependence gives the *quantum state dispersion*. 
To  calculate the  quantum dispersion we  need  to  define the  block *Dispersion*


::

  Dispersion 
  {
    k-path = G-K-M
    number_of_nodes = 10
    k_max = 0.1  
  }

The dispersion of quantum states is calculated at k-points that are nodes of a mesh
in k-space.

The main parameters are:

 ``k-path`` : 
    path of the  direction in  k-space; it  must be  a  string of the  kind G-K-M

 ``number_of_nodes`` : 
   number of nodes in the mesh along each  direction

 ``k-max`` : 
   max value  of  k in  each  direction (default = 0.1)



The dispersion can  be calculated in  general in  a  k-space  dimension between 1 and (3-*simdim*),  where  *simdim*  is  the  simulation  dimension.
Accordingly,  *number_of_nodes*  is    an array  of  dimension equal  to  the k-space  dimension.

For example, for a 2D k-space, a  2D  mesh  is  defined,  by  writing e.g. ::
 
  number_of_nodes = (5,5)


By default, the length in all the directions in   k-space is given by ``k_max``. 
E.g.,  in  a 2D k-space, dispersion is  calculated on a  grid given by  y-axis and z-axis between 0 and ``k_max``

 
If  the  optional keyword ::

  k-path

is  present,  then  the  dispersion is calculated in  1D  along  a  defined  path P1-P2-P3,  for  example G-K-M. Available symmetry points   are G,K,M,M'.
In  this  case,  ``number_of_nodes`` is an  integer.


The dispersion is calculated for each of  the quantum states previously  defined in  the  simulation.








Output
^^^^^^^^^^^^^^^^^^^^^^^

The available output variables for  *Module efaschroedinger* , to be specified in the plot option, are the following:

 ``EigenEnergy`` :
   List of Eigenenergies in eV

 ``ProbabilityDensity``: 
    square modules :math:`|\psi({\bf r})|^2` of the eigenstate wavefunctions 

 ``QuantumDensity``:  
    if  present, quantum  density is  calculated 

 ``k-space_dispersion`` :
    output of  the quantum dispersion

 ``Occupation``:
     probability to find the state occupied. It is calculated assuming Fermi distribution and mean electrochemical potential and temperature:

..  math::
    :nowrap:
    :label:


    \begin{align}
    \bar \mu &= \langle \psi|\mu({\bf r})|\psi \rangle \\
    \bar T &= \langle \psi|T({\bf r})|\psi \rangle 
    \end{align}



If ``ProbabilityDensity`` is specified as plot variable, then ``EigenEnergy`` will plot the levels of the states as constant values on the simulation mesh in addition to the textual file listing all energies.



Module opticskp
----------------------

By defining the **Module** ``opticskp`` , calculation of optical properties is enabled; in particular, 
the optical ``kp`` matrix elements are calculated from the quantum models specified in the **Module**.

The optical spectrum from spontaneous emission is calculated in the following way

..  math::
    :nowrap:
    :label:


    \begin{equation}
         P(\hbar \omega) = \sum_{i,j} \frac{1}{2\pi^2}  \frac{\omega^2_{ij} e^2 }{m^2 c^3}  
         |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
         \frac{\Gamma/2} {(\hbar \omega_{ij} - \hbar \omega)^2 + (\Gamma/2)^2} d\Omega,
    \end{equation}



where :math:`f_i` and :math:`f_j` are the Fermi distributions.

::

  Module opticskp
  {
    name = optics
    regions = quantum
    plot = (optical_spectrum_k_0 )
    initial_state_model = quantum_el
    final_state_model = quantum_hl

    #initial_eigenstates = (0, 9)
    #final_eigenstates = (0, 15)

    polarization = (0, 0, 1)
    Emin = 2.8
    Emax = 3.6
    dE = 0.001
  }


The main parameters are:

 ``initial_state_model``:
     quantum simulations associated to the initial state of the optical transition

 ``final_state_model`` :
    quantum simulations associated to the final  state of optical transition

 ``Emin``:
    minimum  energy 

 ``Emax``:
    max  energy

 ``dE``:
    energy step

 ``plot`` :
    to  select  spectrum in a k-point k0 or  integrated in k-space 


Here, ``initial_state_model`` and ``final_state_model`` are, respectively, the quantum simulations 
( ``efaschroedinger`` module) associated respectively to the initial state of the optical
transition (e.g. electron), and to the final state of optical transition (e.g. hole). 



By  default, all  the eigenstates calculated in  the electron  and  hole quantum simulations are taken in  account  for the  optical calculations.

Optionally, one  can  define  the  range  of  states to  be  used  for  optical  transitions, by  means  of  the  keywords ``initial_eigenstates`` and ``final_eigenstates``. 

A range of energy values must  be  defined in this way::

  Emin = 3.0
  Emax = 5.0
  dE = 0.001

where ``Emin,Emax amd dE`` are  respectively   the  minimum,  maximum and  interval energy. 

Polarization  indicates the  light polarization (vector)..............

By specifying *optical_spectrum_k_0* in  the  ``plot``  statement, the emission  recombination  optical spectrum for  k=0 is  calculated. 


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


The parameters for the *k-space* integration are the following

 ``k_max`` :
    max  value  of k  

 ``quadrature_type`` : 
    type of  integration, default is  *gaussian* 
  

 ``quadrature_order`` : 
    order of  integration, default is  *third*

 ``refine_k_space`` :
    default is *false*, if  *true*, then *adaptive* k-mesh refinement is enabled
 
 ``refine_fraction`` : 
    refinement  parameter
    
 ``relative_accuracy``:
    refinement  tolerance

If ``refine_k_space`` = **true** , that is adaptive k-mesh refinement is enabled, all the elements whose error is greater than the value (1-refine_fraction)* (maximum error) are
going to be refined. In this case, "Error" is just the integrated quantity. The refinement
will end when the tolerance defined by ``relative_accuracy`` is satisfied.


To  run  an  opticskp  sumulation,  one  needs  to  define it in the  solve  statement,  in  this  way::





   Simulation
   { 
      solve = (strain, dd,optics) 
   }


where  *optics*  is  the  name  of  the  defined ``opticskp``  sumulation.
Note  that, in  this  way, the  quantum (efaschroedinger)  simulations needed  for ``opticskp`` are executed but their  results are  not  shown in  the  output.
To plot quantum results you should  explicitly  define   efaschroedinger  simulations in the  solve  statement, in  this  way::




   Simulation
   { 
      solve = (strain, dd, quantum_el,quantum_hl,optics) 
   }





Output
^^^^^^^^^^^^^^^^^^^^^^^


The output variables for optics calculations are:

 ``optical_spectrum_k_0`` : 
    optical emission spectrum for *k=0*.

 ``optical_spectrum`` : 
    optical emission spectrum integrated in  k-space.

 

.. warning::  
            The  two output variables   are  mutually   exclusive.







Example
-----------------------





In this 1D example we will see how to calculate quantum properties of a GaAs/InGaAs **quantum  well**.
:math:`Schr\ddot{o}dinger` equation is solved, with a single-band effective mass model for conduction band and with a 6-band **k.p** model for valence band. Eigenvalues and eigenfunctions are calculated to get energy levels and wavefunctions in the quantum well.
Here is  the  device  definition ::

  Device
  {
   
   meshfile = InGaAs_1D.msh 

   # mesh is drawn in nm
    mesh_units = 1e-9

    Region buffer
    {
      
     material = GaAs

    }


    Region barrier_1
    {
     mesh_regions = (barrier1_cl, barrier1_q)
     
     material =  GaAs
     
    }


    Region QWell
    {
    
     material = InGaAs 
     x = 0.40 

    }


    Region barrier_2
    {
     mesh_regions = (barrier2_q, barrier2_cl)
    
     material = GaAs
    
    }




The InGaAs  well  region (**QWell**) and  the  two  barrier  regions are  collected  in  the  *Cluster*  **Quantum_1**.  Quantum calculation  will  be  restricted  to  this  *Cluster* ::

  Cluster Quantum_1
    {
      regions = (barrier1_q, QWell, barrier2_q)
    } 


First, a strain calculation for the GaAs/InGaAs/GaAs heterostructure, with GaAs as a reference substrate, is performed. ::

  Module elasticity 
   {

    name = strain
     regions = all

     plot = (Strain)


    Physics 
    {
      body_force lattice_mismatch 
      {
        reference_material = GaAs
        structure = zb    
      }
       
    }

    Contact cathode 
    {type = clamp}

  }

See   Module  Elasticity  for  a  detailed  explanation. 

Simulation is performed at equilibrium, so that   we   solve *driftdiffusion* with the option  ::

  coupling = poisson

which  means that only Poisson  equation will  be  solved. ::

  Module driftdiffusion
  { 
   ............
   coupling = poisson
   
   Physics
    {

      strain_simulation = strain

      polarization (piezo, pyro) {}

      recombination srh {}

      recombination direct {}

    }



Finally, we define the model *efaschroedinger*, for quantum effective mass calculations.
We are going to study quantized states of electrons and holes in the quantum well. 
Since the structure is 1D, each eigenstate is characterized by the energy level number *n* and the *k vector* that is perpendicular to the growth direction. In this example, we define two simulations that solve :math:`Schr\ddot{o}dinger` equation for a single *k-vector* (*k = 0*), one for electrons and one for holes. For electrons ::

  Module efaschroedinger

  {
 
   name = quantum_electrons
   regions = Quantum_1 

   plot = (ProbabilityDensity, EigenEnergy)

    
   particle = el
   poisson_simulation = driftdiffusion  #  potential from driftdiffusion 
   strain_simulation = strain 
 

   number_of_eigenstates = 6

   Physics
   {
    model = single_band    
   }

  
  }


Note  that we must  define *poisson* and  *strain*  simulation from  which potential and  strain data have  to  be taken ::

  poisson_simulation = driftdiffusion  #  potential from driftdiffusion 
  strain_simulation = strain 

Electron eigenstates  will  be  obtained from a  single  band model, so we define ::

  Physics
   {
    model = single_band    
   }  



In  a  similar  way,  for  holes ::

  Module efaschroedinger

  {
 
   name = quantum_holes
   regions = Quantum_1 

   plot = (ProbabilityDensity, EigenEnergy)

  
   particle = hl
   poisson_simulation = driftdiffusion  #  potential from driftdiffusion 
   strain_simulation  = strain  
   number_of_eigenstates = 12

   Physics
   {
    model = kp   #  k.p  for  valence  band
    kp_model = 6x6
   }

  
  }

where, in  this  case,  we  use a  6x6 :math:`{\bf k \cdot p}`  model to  calculate  hole  eigenstates ::

  Physics
   {
    model = kp   #  k.p  for  valence  band
    kp_model = 6x6
   }

Finally, in  the  block  *Simulation*,  we  state  the *solve*  order. ::

  solve =  (strain, driftdiffusion, quantum_electrons,quantum_holes ) 


.. warning::  We need to  compute strain and driftdiffusion modules first, to  get  the  correct  parameters  for  :math:`{\bf k \cdot p}`  quantum calculations

Thus, the  block  *Simulation* reads  :: 

  Simulation
  {
   dimension = 1
   temperature = 300
   solve =  (strain, driftdiffusion, quantum_electrons,quantum_holes )
   resultpath =  output  
   output_format = grace
  }

 

Output of  simulation is  shown in the following.

First we look at the output of driftdiffusion calculation.
In Fig. :ref:`Conduction and valence bands<bands>` we have the conduction and valence band profiles in equilibrium condition (contained in the output file *driftdiffusion_msh.dat*); the Fermi level is in correspondance of the zero of energy.


..  _bands :

..  figure:: ../data/bands_fig.png
    :align: center
    :scale: 70%

    Conduction and valence bands




..  _cblevels :

..  figure:: ../data/CB_levels_qw8nm_fig.png
    :align: center
    :scale: 70%

    Conduction band levels


..  _vblevels :

..  figure:: ../data/VB_levels_qw8nm_fig.png
    :align: center
    :scale: 70%

    Valence band levels


Then in Fig. :ref:`Conduction band levels<cblevels>` we  show
the conduction band profile and the the quantized levels for electrons (from the file *quantum_electrons_msh.dat*) in the InGaAs quantum well; for this 8 nm-wide quantum well only the first two energy levels are confined. 


The first 12 quantized levels for the holes (from the file *quantum_holes_msh.dat*) in the InGaAs quantum well are shown in Fig. :ref:`Valence band levels<vblevels>` (there are 6 couples of degenerated energy levels). 






Fig. :ref:`CB wavefunctions<cbwave>` show the wavefunctions (square module) for the energy states in the conduction band: only the first two states are confined in the conduction band and are shown here.

 
Finally, in Fig. :ref:`VB wavefunctions<vbwave>`
there are the wavefunctions for the first 10 confined energy states in the valence band. The first of each couple of degenerated states is shown.
States from the heavy hole and the light hole bands are visible. 







..  _cbwave :

..  figure:: ../data/efa_ex1_wave_el.png
    :align: center
    :scale: 70%

    CB wavefunctions   



..  _vbwave :

..  figure:: ../data/wave_hl_ex1.png
    :align: center
    :scale: 70%

    VB wavefunctions   






.. rubric:: Footnotes


..   </marker>


.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%



