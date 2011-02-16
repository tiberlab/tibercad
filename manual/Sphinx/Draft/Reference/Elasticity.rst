.. _Elasticity:

Elasticity
=================================================


Stiffness/Isotropic
Kind: bulk 
Name: Stiffness
Type: isotropic
database_section: none
Options: 
young  (double,0.0)
poisson (double [-0.5,1.0],0.0)

Stiffness/Anisotropic
Kind: bulk 
Name: Stiffness
Type: anisotropic
database_section = stiffness
Options:
none

Body Force/constant
Kind: bulk 
Name: BodyForce
Type: constant
database_section: none
Options: 
force (double vector,(0.0,0.0,0.0))

Body Force/Lattice Mismatch
Kind: bulk 
Name: BodyForce
Type: lattice_mismatch
database_section: none
Options:
x (double [0.0,1,0],0.0)
x-growth-direction (double vector)
y-growth-direction (double vector)
z-growth-direction (double vector)

Body Force/Converse
Kind: bulk 
Name: BodyForce
Type: converse
database_section: piezoelectricity
Options: 
poisson_simulation (string,”none“)

Boundary/clamp
Kind: interface
Name: Contact
Type: converse
database_section: none
Options: 
poisson_simulation (string,”none“)

Boundary/Surface Force
Kind: interface
Name: Contact
Type: surface_force
database_section: none
Options: 
force (double vector,(0.0,0.0,0.0)

