.. _Glossary:

Glossary
========

..

.. glossary::
   :sorted:
   
   Stiffness;Isotropic

      | Level : 
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

      | Level : 
      | Kind: bulk 
      | Name: Stiffness
      | Type: anisotropic
      | database_section = stiffness
      | Options:
      |     none
      | 
      | Description :

   Body Force;constant

      | Level : 
      | Kind: bulk 
      | Name: BodyForce
      | Type: constant
      | database_section: none
      | Options: 
      |      force (double vector,(0.0,0.0,0.0))
      | 
      | Description :


   Body Force;Lattice Mismatch

      | Level : 
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

      | Level : 
      | Kind : bulk 
      | Name : BodyForce
      | Type : converse
      | database_section : piezoelectricity
      | Options : 
      |     poisson_simulation (string,"none")
      | 
      | Description :


   Boundary;clamp

      | Level : 
      | Kind: interface
      | Name: Contact
      | Type: converse
      | database_section: none
      | Options: 
      |     poisson_simulation (string,"none")
      | 
      | Description :


   Boundary;Surface Force

      | Level : 
      | Kind: interface
      | Name: Contact
      | Type: surface_force
      | database_section: none
      | Options: 
      |     force (double vector,(0.0,0.0,0.0)
      | 
      | Description :

