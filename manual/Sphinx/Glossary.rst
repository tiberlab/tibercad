.. _Glossary:

Glossary
========

.. glossary::
   :sorted:
   
   Stiffness;Isotropic
      | Level:
      | Kind: bulk 
      | Name: Stiffness
      | Type: isotropic
      | database_section: none
      | Options: 
      |      young  (double,0.0)
      |      poisson (double [-0.5,1.0],0.0)
      | 
      | Description :
      | 


   Stiffness;Anisotropic
      | Level:
      | Kind: bulk 
      | Name: Stiffness
      | Type: anisotropic
      | database_section = stiffness
      | Options:
      |     none
      | 
      | Description :
      | 
      
   Body Force;constant
      | Level:
      | Kind: bulk 
      | Name: BodyForce
      | Type: constant
      | database_section: none
      | Options: 
      |      force (double vector,(0.0,0.0,0.0))
      | 
      | Description :
      | 

   Body Force;Lattice Mismatch
      | Level:
      | Kind: bulk 
      | Name: BodyForce
      | Type: lattice_mismatch
      | database_section: none
      | Options:
      |     x (double [0.0,1,0],0.0)
      |     x-growth-direction (double vector)
      |     y-growth-direction (double vector)
      |     z-growth-direction (double vector)
      | 
      | Description :
      | 

   Body Force;Converse
      | Level:
      | Kind: bulk 
      | Name: BodyForce
      | Type: converse
      | database_section: piezoelectricity
      | Options: 
      |     poisson_simulation (string,"none")
      | 
      | Description :
      | 

   Boundary;clamp
      | Level:
      | Kind: interface
      | Name: Contact
      | Type: converse
      | database_section: none
      | Options: 
      |     poisson_simulation (string,"none")
      | 
      | Description :
      | 

   Boundary;Surface Force
      | Level:
      | Kind: interface
      | Name: Contact
      | Type: surface_force
      | database_section: none
      | Options: 
      |     force (double vector,(0.0,0.0,0.0)
      | 
      | Description :
      | 

   Solvers;nonlinear
      | Level: 
      | Kind:
      | Name: 
      | Type: solvers
      | database_section: none
      | Options: 
      |     petsc (string) : uses the PETSc nonlinear solver (SNES) (default)
      |               relative_tolerance : convergence criterion based on relative residual l_2-norm
      |               absolute_tolerance : convergence criterion based on the l_2-norm of the residual
      |               max_iterations : maximum number of iterations
      |               step_tolerance : tolerance criterion based on the l_2-norm of the correction step
      |               max_step : maximum linesearch step (l_2-norm)
      |               divergence_tolerance : divergence criterion
      |               
      |     linesearch (string) : uses a linear linesearch implemented in TiberCAD
      |               absolute_tolerance : convergence criterion based on the l_2-norm of the residual
      |               step_tolerance : tolerance criterion based on the l_infinity-norm of the correction step
      |               max_iterations : maximum number of iterations
      | 
      | Description : Nonlinear solvers  are based on iterative methods, solving in each iteration a linear system.
      | 
      
   Solvers;linear
      | Level: 
      | Kind:
      | Name: 
      | Type: solvers
      | database_section: none
      | Options: 
      |     petsc (string) : uses the PETSc linear solver (SNES) (default)
      |               method : the Krylov subspace method to be used     
      |                  bcgs =  BiCGSTAB
      |                  bcgsl =
      |                  gmres =  Generalized Minimal Residual
      |                  bicg =  BiConjugate Gradient
      |                  cg =  Conjugate Gradient
      |                  cgs =  Conjugate Gradient Squared
      |                  richardson =  Richardson
      |                  pconly = only apply preconditioner
      |     
      |     relative_tolerance : convergence criterion based on relative residual l_2-norm
      |     
      |     absolute_tolerance : convergence criterion based on the l_2-norm of the residual
      |     
      |     max_iterations : maximum number of iterations
      |     
      |     preconditioner (string) : attribute of pconly
      |               method : the Krylov subspace method to be used     
      |                  lu =  LU
      |                  ilu = incomplete LU
      |                  jacobi =  Jacobi
      |                  cholesky =  Cholesky
      |                  none =  no preconditioning
      |                  composite =  use a combination of Jacobi and iLU
      | 
      | Description : Linear solvers  are based on iterative methods, solving in each iteration a linear system.
      |  
      
   Simulations;sweep
      | Level: Sub-Block
      | Kind:
      | Name: sweep
      | Type: simulation block
      | database_section: none
      | Options: 
      |   solve = a list of simulations to be solved, eg. (strain, dd)
      |   variable =  the sweep variable, eg. $Vd
      |   start =  the first value of the sweep
      |   stop =  the last value of the sweep
      |   steps =  the number of steps the interval (start, stop) gets subdivided in
      |   values =  instead of start, stop and steps, provide the sweep values explicitly
      |   plot_data =  set to true if you want to plot after each step the simulation results (default is false)
      |   min_step =  the minimum absolute step size
      |   max_step = the maximum absolute step size
      |   initial_step = the absolute initial step size
      |   min_relative_step =  minimum relative step size, default is 1e-3
      |   max_relative_step =  maximum relative step size, default is 1
      |   initial_relative_step =  initial relative step size, default is 1
      | 
      | Description : Performs a linear sweep for a given variable.
      | 
      
   Solvers;selfconsistent
      | Level: Sub-Block
      | Kind:
      | Name: selfconsistent
      | Type: solvers
      | database_section: none
      | Options: 
      |   max_iteration =  the maximum number of iterations. Currently, when the maximum number of iterations is reached, the program only issues a warning and proceeds.
      |   relative_tolerance =  the relative convergence tolerance for the observed variable in terms of the l2-norm
      |   absolute_tolerance =  the absolute convergence tolerance for the observed variable in terms of the maximum-norm
      |   relaxation_factor =  an optional relaxation factor (default is 1) to be applied to the observed variable
      |   solve =  the simulations to be solved
      | 
      | Description : Solves models in an iterative way to obtain a selfconsistent solution, optionally using a relaxation approach.
      | 

   Heat Source;constant
      | Level: Sub-Block
      | Heat Source/constant
      | Kind: bulk 
      | Name: HeatSource
      | Type: constant
      | database_section: none
      | Options: 
      |   H (double,0.0)
      | 
      | Description : 
      | 

   Heat Source;joule
      | Level: Sub-Block
      | Heat Source/joule
      | Kind: bulk 
      | Name: HeatSource
      | Type: joule
      | database_section: none
      | Options: 
      |   dd_simulation (string,"driftdiffusion")
      | 
      | Description : 
      | 

   Boundary;Heat reservoir
      | Level: Sub-Block
      | Boundary/Heat reservoir
      | Kind: interface
      | Name: Contact
      | Type: heat_reservoir
      | database_section: none
      | Options: 
      |   temperature (positive double,300)
      | 
      | Description : 
      | 

   Boundary;Surface Thermal Resistance
      | Level: Sub-Block
      | Boundary/Surface Thermal Resistance
      | Kind: interface
      | Name: Contact
      | Type: surface_resistance
      | database_section: none
      | Options: 
      |   r_surf (positive double 0.0)
      |   temperature (positive double (K),300)
      | 
      | Description : 
      | 

   Boundary;Flux
      | Level: Sub-Block
      | Boundary/Flux
      | Kind: interface
      | Name: Contact
      | Type: flux
      | database_section: none
      | Options: 
      |   Flux  (double vector (W/cm2, (0.0,0.0,0.0))
      | 
      | Description : 
      | 
      
   Module efaschroedinger
      | The EFA calculation  of eigenstates and  eigenfunctions are performed by the **Module** efaschroedinger.
      | 
      | Options: 
      |   particle = el | ..
      |      particle = string defines for which particle (electron or hole) Schroedinger equation is solved Possible values are el and hl. A different Module efaschroedinger has to be defined for each particle to be solved.
      |
      |   poisson_model_name = dd | ..
      |      poisson_model_name = string defines the name of the simulation (Module driftdiffusion) that can provide electric potential
      |
      |   strain_model_name = strain | macrostrain
      |      strain_model_name = string defines the name of the simulation (Module macrostrain) that can provide elastic strain
      |
      |   name = quantum_el | ..
      |
      |   regions = quantum | ..
      |      regions = string defines the regions associated to this EFA simulation
      |
      |   plot = EigenFunctions | EigenEnergy | EnergyLevels
      | 
      | can include a Solver Sub-Block 
      | 
      |    Options: 
      |       number_of_eigenstates = 10 | 30 
      | 
      | can include a Physics Sub-Block 
      | 
      |    Options: 
      |       model = conduction_band | ..
      | 

      
   Module quantumdispersion
      | With the Module quantumdispersion it is possible to calculate the dependence of quantum eigenstates on **k** -vector. Such dependence gives the *quantum state dispersion* . The simulation name is **quantumdispersion** .
      | 
      | Options: 
      |   simulation_name = dispersion1D_el | ..
      |      
      |
      |   regions = all | ..
      |      
      |
      |   quantum_simulation = quantum_el | ..
      |      quantum simulation is the name of the efaschroedinger simulation.
      |
      |   min_eigenvalue_number = [int] i.e. 0
      |      the dispersion is calculated for the states number **i** , where i is included between the values of min and max eigen numbers.
      |
      |   max_eigenvalue_number = [int] larger than min_eigenvalue_number i.e. 5
      |      the dispersion is calculated for the states number **i** , where i is included between the values of min and max eigen numbers.
      |
      |   wedge = half
      |      defines properties of K-space.
      |
      |   k_space_dimension = 1
      |      defines properties of K-space.
      |
      |   k1 = (0, 0.1, 0)
      |      defines properties of K-space.
      |
      |   number_of_nodes = (10)
      |
      |   output_format = grace
      |      the main format to view data with Grace software
      |
      |   plot = k-space_dispersion | ..
      |      variables to include in the plot graphs
      | 

      
   Module quantumdensity
      | In Module quantumdensity you can define the calculation of particle (electron,hole) density, 
      | based on the result of a previous calculation of the system eigenstates (e.g. with
      | efaschroedinger module). Quantum density may be obtained with an analytical or a
      | numerical calculation.
      | 
      | Options: 
      |   name = dens_el | ..
      |      
      |
      |   regions = quantum | ..
      |      
      |
      |   quantum_simulation = quantum_el | ..
      |      quantum simulation is the name of the efaschroedinger simulation.
      |
      |   initial_eigenstates_number = [int] i.e. 10
      |      initial number of eigenstates for the Schroedinger equation
      |
      |   degeneracy = [int] i.e. 2 [int] 
      |      degeneracy of the quantum state
      |
      |   refine_fraction = [double] i.e. 0.20
      |      
      |   relative_accuracy = [double] i.e. 0.01
      |      
      |   refine_k_space = [bool] i.e. true
      |      
      |   output_density_in_k_space = [bool] i.e. true
      |      
      |   uniform_refinement = [bool] i.e. false
      |      
      |   mesh_order = FIRST | ..
      |      
      |   first_state = [int] i.e. 3
      |
      |   analytic = [bool] i.e. true
      |      If true then the density is calculated analytically, otherwise numerically.
      |
      |   wedge = half
      |      defines properties of K-space.
      |
      |   k_space_dimension = 0
      |      defines properties of K-space.
      |
      |   k1 = (0, 0.1, 0)
      |      defines properties of K-space.
      |
      |   number_of_nodes = (4)
      |
      |
      |   plot = quantum_density | ..
      |      variables to include in the plot graphs
      | 
      
   Module opticskp
      | By defining the Module **opticskp** , calculation of optical properties is enabled; in particular, 
      | the optical kp matrix elements are calculated from the quantum models specified in the Module.
      | 
      | Options: 
      |   name = optics
      |      
      |
      |   regions = quantum
      |      
      |
      |   initial_state_model = quantum_el 
      |      quantum simulation associated to the initial state of optical transition (e.g. electron)
      |
      |   final_state_model = quantum_hl
      |      quantum simulation associated to the final state of optical transition (e.g. hole)
      |
      |   initial_eigenstates = (0, 9)
      |      referred to the range of eigenstates to be taken in account for optical calculations
      |
      |   final_eigenstates = (0, 15)
      |      referred to the range of eigenstates to be taken in account for optical calculations
      |
      |   polarization = (0, 0, 1)
      |      
      |
      |   Emin = [int] i.e. 2.8
      |      minimum value for the range of energy
      |
      |   Emax = [int] greater than Emin i.e.3.6
      |      maximum value for the range of energy
      |
      |   dE = 0.001
      |      incresing step selected for the range of energy
      |
      |   plot = (optical_spectrum_k_0 ) | ..
      |      optical emission spectrum for *k=0*
      | 
      
   Module opticalspectrum
      | By defining the Module **opticskp** , calculation of optical properties is enabled; in particular, 
      | the optical kp matrix elements are calculated from the quantum models specified in the Module.
      | 
      | Options: 
      |   k_space_dimension = 2
      |        k_space_dimension = **1** for 2D simulations, **2** for 1D simulations. k-space basis is
      |        **true** if the k-space is defined by means of k-vectors; if **false** , vectors are expressed in
      |        real space
      |
      |   k1 = (0, 0, 0.1)
      |   k2 = (0, 0.1, 0)
      |
      |   refine_fraction = 0.30 
      |      
      |
      |   relative_accuracy = 0.01
      |      
      |
      |   refine_k_space = true
      |      **true** , that is adaptive k-mesh refinement is enabled, all the elements 
      |      whose error is greater than the value (1-refine fraction) (maximum error) are
      |      going to be refined. In this case, Error is just the integrated quantity. 
      |      The refinement will end when the *relative_accuracy* is obtained.
      |
      |   number_of_nodes = (2, 2)
      |      number of elements in k mesh, along each direction
      |
      |   wedge = half | quarter
      |      to reduce calculation time, by exploiting symmetry
      |   
      |   optical_matr_elem_model = opticskp
      |      name of the *opticskp* model associated
      |   
      |   polarization = (0, 0, 1)
      |      light polarization (vector)
      |   
      |   Emin = [int] i.e. 3.0
      |      minimum value for the range of energy
      |
      |   Emax = [int] greater than Emin i.e. 5.0
      |      maximum value for the range of energy
      |
      |   dE = 0.001
      |      incresing step selected for the range of energy
      |
      |   plot = (optical_spectrum  | .. )
      |      k-space integrated optical emission spectrum
      | 
      
