

Elasticity
=================================================


Introduction
----------------------

Elasticity is a Finite Elements solver for mechanical equilibrium problems. 
It brings features typically developed for structural stability solvers into device modeling. 
The coupled treatment of the electro-mechanical problem within a unique framework results very useful to explore for multidisciplinary ideas. 


Details about the theory of continuous elasticity applied to device mechanical deformation 
can be found in [Povolotskiy]_ .

TiberCAD computes the mechanical deformation of a body subjected to external forces by means of the equilibrium mechanical equation, i.e.

.. math::
   :label: ela_eq_001


    \frac{\partial \sigma_{ij}}{\partial x_j} = f_i


wheres :math:`\sigma` is the stress second-rank tensor and f is the external body force applied to the system. 
As we will see below, any strain and stress source can be appropriately mapped into an equivalent body force. 

The stress tensor :math:`\sigma` is related to the mechanical strain by means of the constitutive relationships :math:`\sigma_{ij}=C_{ijlk}\epsilon_{lk}` 
where **C** is the stiffness fourth-rank tensor. 
The strain is related to the displacement u by the expression

.. math::
   :label: ela_eq_002

    
    \epsilon_{kl}=\frac{1}{2}\left(\frac{\partial u_l}{x_k}-\frac{\partial u_k}{\partial u_l} \right )
    

Relying on the symmetry of the strain and stiffness tensor the final equation reads as

.. math::
   :label:

    
    \frac{\partial}{\partial x_j}C_{ijlk}\frac{\partial u_l}{\partial x_k} = f_i
    

The non-linear strain is computed by applying the mechanical equilibrium equation on the deformed mesh.  
The number of the shape iteration is set by the keyword ``shape_iteration`` (default = 0, i.e. no deformation is performed) 
whereas the maximum tolerated error can be indicated with shape_error ( default = :math:`1e^{-3}` ).




.. _ElasticityGetting:

..  GETTING STARTED PART


Example
----------------------

Assuming a small displacement, Elasticity computed the strain by solving the equation

.. math::
   :label: ela_eq_011
   
   
   \frac{\partial}{\partial x_j}C_{ijlk}\frac{\partial u_l}{\partial x_k} = f_i
   

where ``C`` is the stiffness constant  and ``f`` is the total external body force. 
The latter may include several contribution such as the converse piezoelectric effect the thermal expansion and a lattice mismatch induced strain. 
The latter, described in the example below, occurs when an interface is created between two materials with different lattice constants. 
Let :math:`\epsilon^{LM}` be the lattice mismatch strain, then the effective body force reads as :  

.. math::
   :label: ela_eq_012

    
    f_i = -\frac{\partial }{\partial x_j} C_{ijlk}\epsilon_{lk}^{LM}
    



In the following we will compute the strain induced by the lattice mismatch in a system comprising a layer of GaN between two contacts of :math:`Al_{x}Ga_{1-x}N` . 
Once the mesh is created with ``gmsh`` (distributed along TiberCAD ) we may start to build the input file. 
First, we have to insert the region section ::

  Device
  {
    dimension = 3 
    meshfile = mesh.msh
    mesh_units = 1e-9
    material = AlGaN
    x = 0.2
    x-growth-direction = (-1,0,1,0) 
    y-growth-direction = (-1,2,-1,0) 
    z-growth-direction = (0,0,0,1)
    Region Well  {material = GaN}
  }

All the options included in this section, such ad the material, axis growth and alloy concentration, apply for the whole structure.  
If we want to specify a region with different properties, we may use the keyword ``Region`` and refer to a region name, as specified in the mesh file. 

The second part of the input file is devoted to the module declaration ::

  Module elasticity 
  {
    Physics 
   {
     BodyForcelattice_mismatch 
     {
       x-growth-direction = (-1,0,1,0) 
       y-growth-direction = (-1,2,-1,0) 
       z-growth-direction = (0,0,0,1)
       reference_material = AlGaN
       x = 0.2 
      }
     Contact Base {type = clamp}
     }
    }

In ``physics`` we declare the physical models to be applied to our device. 
In our case, with ``BodyForce`` lattice mismatch we are adding a body force into the system, 
induced by the lattice mismatch with the reference lattice which in this case is :math:`Al_{0.2}Ga_{0.8}N` . 
In order to avoid a free standing device we may want to freeze a surface. With ``Boundary Clamp`` we fix, in this case, the region Base.

The third part is simply Simulation ``{solve = elasticity}`` where the default name *elasticity* is used. 
By default, files for 3D system are written in vtu which can be read from Paraview. 
Output data about strain are shown below.

.. _fig.elasticity01: 

.. figure:: ../data/elasticity01.png
   :align: center
   :scale: 80%
   
   The 3D output

   
   
.. _fig.elasticity02: 

.. figure:: ../data/elasticity02.png
   :align: center
   :scale: 80%
   
   Strain effect on directions


..   EOF GETTING STARTED PART



.. rubric:: Footnotes




