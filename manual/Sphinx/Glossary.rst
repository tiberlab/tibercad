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

