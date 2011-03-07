..   <marker>

.. _ElasticityReferenceguide:

Elasticity
=================================================

.. index:: double:Stiffness;Isotropic

| Stiffness/Isotropic
| Kind: bulk 
| Name: Stiffness
| Type: isotropic
| database_section: none
| Options: 
| young  (double,0.0)
| poisson (double [-0.5,1.0],0.0)

.. index:: double:Stiffness;Anisotropic

| Stiffness/Anisotropic
| Kind: bulk 
| Name: Stiffness
| Type: anisotropic
| database_section = stiffness
| Options:
| none

.. index:: double:Body Force;constant

| Body Force/constant
| Kind: bulk 
| Name: BodyForce
| Type: constant
| database_section: none
| Options: 
| force (double vector,(0.0,0.0,0.0))

.. index:: double:Body Force;Lattice Mismatch

| Body Force/Lattice Mismatch
| Kind: bulk 
| Name: BodyForce
| Type: lattice_mismatch
| database_section: none
| Options:
| x (double [0.0,1,0],0.0)
| x-growth-direction (double vector)
| y-growth-direction (double vector)
| z-growth-direction (double vector)

.. index:: double:Body Force;Converse

| Body Force/Converse
| Kind: bulk 
| Name: BodyForce
| Type: converse
| database_section: piezoelectricity
| Options: 
| poisson_simulation (string,"none")

.. index:: double:Boundary;clamp

| Boundary/clamp
| Kind: interface
| Name: Contact
| Type: converse
| database_section: none
| Options: 
| poisson_simulation (string,"none")

.. index:: double:Boundary;Surface Force

| Boundary/Surface Force
| Kind: interface
| Name: Contact
| Type: surface_force
| database_section: none
| Options: 
| force (double vector,(0.0,0.0,0.0)


..   </marker>
