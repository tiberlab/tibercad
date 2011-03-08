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

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

      
   Parameter
   
      | Level: 
      | Kind:
      | Name: 
      | Type: 
      | database_section: none
      | Options: 
      |     
      | 
      | Description :

