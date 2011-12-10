..   <marker>

..  _ElasticityTheory:

Elasticity
=================================================




Theory
-------------


Elasticity is a Finite Element solver for mechanical equilibrium problems. 
It brings features typically developed for force balance problems into device modeling. 
The coupled treatment of the electro-mechanical problem within a unique framework results very useful to explore for multidisciplinary ideas. 
Details about the theory of continuous elasticity applied to device mechanical deformation 
can be found in [Povolotskyi]_ .

TiberCAD computes the mechanical deformation of a body subjected to external forces by means of the equilibrium mechanical equation, i.e.

.. math::
   :label: mec_eq

    
    \frac{\partial \sigma_{ij}}{\partial x_j} = f_i
    

wheres :math:`\sigma` is the stress second-rank tensor and f is the external body force applied to the system. 
As we will see below, any strain and stress source can be appropriately mapped into an equivalent body force. 

The stress tensor :math:`\sigma` is related to the mechanical strain by means of the constitutive relationships :math:`\sigma_{ij}=C_{ijlk}\epsilon_{lk}` 
where **C** is the stiffness fourth-rank tensor. 
The strain is related to the displacement u by the expression

.. math::
   :label: strain_eq 

    
    \epsilon_{kl}=\frac{1}{2}\left(\frac{\partial u_l}{x_k}-\frac{\partial u_k}{\partial u_l} \right )
    

Because of the symmetry of the strain and stiffness tensor the final equation reads as

.. math::
   :label: strain_eq2

    
    \frac{\partial}{\partial x_j}C_{ijlk}\frac{\partial u_l}{\partial x_k} = f_i
    

the computed displacements are used to deform the mesh. 



Stiffness constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default the stiffness constant is anisotropic. For the zincoblend crystal structure the independent coefficients are 
(with the Voigt notation) :math:`C_{11}, C_{12}, C_{44}`. 

..  math::
    :label: stiffness
    
   
    C = \left(
    \begin{array}{cccccc}
    C_{11} & C_{12} & C_{12} & 0 & 0 & 0\\
    C_{12} & C_{11} & C_{12} & 0 & 0 & 0\\
    C_{12} & C_{12} & C_{11} & 0 & 0 & 0\\
    0 & 0 & 0 & C_{44} & 0 & 0\\
    0 & 0 & 0& 0 & C_{44} & 0\\
    0 &  0  &  0 & 0 & 0 & c_{44} 
    \end{array}
    \right)
   
    
|

For the wurtzite structure we have :math:`C_{11}, C_{12}, C_{13}, C_{44}`.

..  math::
    :label: stiffness2
    
   
    C = \left(
    \begin{array}{cccccc}
    C_{11} & C_{12} & C_{13} & 0 & 0 & 0\\
    C_{12} & C_{11} & C_{12} & 0 & 0 & 0\\
    C_{13} & C_{12} & C_{11} & 0 & 0 & 0\\
    0 & 0 & 0 & C_{44} & 0 & 0\\
    0 & 0 & 0& 0 & C_{44} & 0\\
    0 &  0  &  0 & 0 & 0 & \frac{C_{11} - C_{22}}{2}   
    \end{array}
    \right)
   

|

While the anisotropic model is included by default, the ``isotropic`` one must be explicitly
indicated.

An isotropic stiffness can be defined within the block ``Stiffness`` with the keyword ``isotropic``. 
In this case, the only independent parameters are the Young modulus :math:`E`
(``young`` in GPa) and the Poisson's ratio :math:`\nu` (``poisson``). 

..  math::
    :label: anisotrop_stiff
    
   
    C =
    \frac{E}{(1+\nu)(1-2\nu)} \left(
    \begin{array}{cccccc}
    1-\nu & \nu & \nu & 0 & 0 & 0\\
    \nu & 1-\nu & \nu & 0 & 0 & 0\\
    \nu & \nu & 1-\nu & 0 & 0 & 0\\
    0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0\\
    0 & 0 & 0& 0 & \frac{1-2\nu}{2} & 0\\
    0 &  0  &  0 & 0 & 0 & \frac{1-2\nu}{2}   
    \end{array}
    \right)
   

|


Example::

  Stiffness isotropic
    {
     young = 129
     poisson = 0.349
    }




Module options
---------------------

The following options influence the behaviour of the Elasticity module:

 ``non_linear_strain`` : bool
   if :math:`true` the strain is computed iteratively until the convergence on the structure deformation is reached  (default = :math:`false`) 

 ``shape_iterations`` : integer
    defines the number of  shape iterations for non-linear strain computation. 
    The default is ``0``, meaning that the force balance problem is solved only once. 

 ``shape_error`` : double
   defines the maximum tolerated error in non-linear strain computation (default = :math:`1e-3`) 
  

Solution/Plot variables
-----------------------

The solution variables available for plotting and for interaction with other modules are
given in :ref:`Plotting variables<el_solutions>`.




Solver section
--------------------

The ``Solver`` section of the ``Elasticity``   module refers to a linear solver.
See section :ref:`Linear_solver`  for details on linear solver options.

 

Physics section
--------------------

In the following we describe all the submodels. As mentioned in the Introduction  in section :ref:`InputFileGetting`, submodels can be restricted to a subset of simulation regions.

The user can specify the following physical  models:

  * ``Body  force``

  * ``Stiffness``

Body  force
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The  Body   force  model   implements  an  external  force  applied  to  the  device.


The  keyword  is ``body_force`` , i.e ::

  body_force [type]
  {
  }

The type can be ``constant``, ``converse_piezo``   or ``lattice_mismatch``.


Constant
...............


A constant value body force can be included by means of the keyword ``constant`` .

Example::

  body_force constant 
    {
     F = 1
    }


F  is  a  force in :math:`N/m^3`    


Lattice mismatch
.................

When two crystals with different crystal structure are put in contact the lattice mismatch between them may induce a strain  :math:`\epsilon^{LM}`  and, therefore, a stress :math:`\sigma =C\epsilon^{LM}` . 

This contribution can be mapped in the body force

.. math::
   :nowrap:
   :label:

    \[
    f_i = -\frac{\partial }{\partial x_j} C_{ijlk}\epsilon_{lk}^{LM}
    \]


and  can  be included with the keyword ``lattice_mismatch`` .

The strain source can be computed only once the reference lattice is identified. 
The reference material and its growth axis must be included following the same syntax of the device section (see :ref:`InputFileGetting`).

Example::

   body_force lattice_mismatch
    {
     reference_material = AlGaN
     structure = wz
     x = 0.2
     x-growth-direction = (1,0,0,0)
     y-growth-direction = (0,1,-1,0)
     z-growth-direction = (0,0,0,1)
    }




.. note::
          
        When a lattice mismatch force is included in the calculation, the displacement will be computed with respect the reference lattice.


Converse  piezo
.................


In presence of an electric field there might develop an additional stress source due to the so-called converse piezoelectric effect, given by :

.. math::
   :nowrap:
   :label:

    \[
    \sigma_{il}^{SC}=-e_{ijk}E_k
    \]

The converse piezoelectric effect can be included with the keyword ``converse_piezo`` and an electrostatic simulation must be indicated. 

The effective body force reads as :

.. math::
   :nowrap:
   :label:

    \[
    f_i = -\frac{\partial}{\partial x_j} \sigma_{ij}^{CP}
    \]


Example::

   body_force  converse_piezo  
    {
     DD_simulation = dd
    }


where  ``DD_simulation``  indicates the relevant electrostatic simulation.



Stiffness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



The  keyword  is ``stiffness`` , i.e ::

  stiffness [type]
  {
  }


The type can be ``isotropic``   or ``anisotropic``.  The  latter is  the default.

isotropic
.................

When the stiffness constant is  ``isotropic``, the  user  has  to  provide the  Young module (``young``) and Poisson ratio (``poisson``).

Example::

   stiffness  isotropic 
    {
        young = 129
        poisson = 0.349
    }



anisotropic
.................

Even  when the  ``anisotropic`` constant is  chosen, it  is  possible to override it in the input file.

Example::

   stiffness  anisotropic 
    {
        C11 =55
        C22 = 123
        C44 = 103
         
    }





Boundary conditions
-------------------

Different boundary conditions can be implemented for Elasticity simulations, by  means  of  the  block  ``Contact``.

  * ``Surface force``

  * ``Clamp``

  * ``Custom``


Surface force
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Surfaces forces :math:`f^0` are applied by imposing :math:`\sigma_{ij} n_{j} = f_i^0` along the surface with normal :math:`n`. 
This boundary condition can be used with the keyword ``surface_force``.  
An example is shown below. ::

  Contact base
    {
     type = surface_force
     force = (0,0,0.5) 
    }

where ``force`` is  the applied force in  GPa.


Clamp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



On the other hand, one may want to fix some surface of the device. The keyword ``clamp`` freezes all nodes of a given surface. 

Example::

  Contact substrate 
    {
     type = clamp
    }







Custom
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``Custom`` imposes a general  boundary  condition by means of the matrix and vector constrain H and R respectively, which are  related to the displacement by means of  :math:`Hu=R`.

Example::

 
   Contact substrate 
     {

      type = Custom
      H = ((1,0,0),(0,1,0),(0,0,1))
      R = (1,0,1)

     }



..  _el_solutions :

..  math::
    :nowrap:
    :label:
    
     \begin{table}[!ht]
     \center
     \begin{tabular}{l|c|l}
     \multicolumn{3}{c}{\textbf{Solution Table}} \\
     \hline
     \textbf{Keyword}  & \textbf{Description} & \textbf{Units}  \\
     \hline
     \hline
     \texttt{Strain} & total strain  & - \\
     \texttt{StrainCrystal} & total strain in the crystal system & - \\
     \texttt{Stress} & total stress & GPa \\
     \texttt{StressCrystal} & total stress in the crystal system & GPa \\
     \texttt{Displacement} & displacement  & m  \\
     \texttt{ForceSource} & force   & N/m3  \\
     \texttt{StrainSource} & strain source  & -  \\
     \texttt{StressSource} & stress source   & GPa \\
     \end{tabular}
     \caption{Solution/Plot variables}
     \label{table:el_solutions}
     \end{table}

|



Example 1: Strain  in  heterostructures
-------

In the following example we will compute the strain induced by the lattice mismatch in a system comprising a layer of GaN between two contacts of :math:`Al_{x}Ga_{1-x}N` . 
Once the mesh is created with *gmsh* (distributed along with TiberCAD ) we may start to build the input file. 
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
    Region Well
      {material = GaN}

    Region Barr1
      {
      }

    Region Barr2
      {
      }


  }

All the options included in this section, such as the material, axis growth and alloy concentration, apply for all  the  Regions defined in the whole structure.  
If we want to specify  different properties for  a  Region, we may define them in  the relevant Region block. E.g. ::
 
  Region Well
      {material = GaN}

 



The second part of the input file is devoted to the module declaration ::

  Module elasticity 
  {
    Physics 
    {
      body_force lattice_mismatch 
      {
        x-growth-direction = (-1,0,1,0) 
        y-growth-direction = (-1,2,-1,0) 
        z-growth-direction = (0,0,0,1)
        reference_material = AlGaN
        x = 0.2 
       }
    
     }

    Contact Base 
        {type = clamp}

   }


In  ``Physics`` we declare the physical models to be applied to our device. 
In our case, with ``body_force lattice_mismatch`` we are adding a body force into the system, 
induced by the lattice mismatch with the reference lattice which in this case is :math:`Al_{0.2}Ga_{0.8}N` . 
In order to avoid a free standing device we may want to freeze a surface. With the Boundary type ``Clamp`` we fix, in this case, the region Base.


.. _fig.elasticity01: 

.. figure:: ../data/elasticity01.png
   :align: center
   :scale: 100%
   
   The 3D output

   
   
.. _fig.elasticity02: 

.. figure:: ../data/elasticity02.png
   :align: center
   :scale: 100%
   
   Strain effect on directions


The third part is simply Simulation ``{solve = elasticity}`` where the default name ``elasticity`` is used. 
By default, files for 3D system are written in *vtu* format which can be read from Paraview. 
Output data about strain are shown below.






Example 2: Piezoelectric nanogenerator
------------------------------------

In this example we will compute the output potential of a piezoelectric nanogenerator based on a vertical compressed ZnO nanowire [1,2]. The cylindrical columni, oriented along the :math:`z`-axis, has a radius of 150 nm and is 4 :math:`\mu m` high. 
In the device section we specify the name of the region (*Column*) and the doping. ::

  Device
   {

    meshfile = mesh.msh
    mesh_units = 1e-6

    material = ZnO

    dimension = 3

        Region Column
        {

          material = ZnO

          Doping
              {
                Nd = 1e16
                type = donor
                Ed = 0.035
              }

         }

  }
  

In the Physics section we use define isotropic elastic constant by means of the Young module and Poisson ratio. ::

  Module elasticity
  {

   plot = (Strain,Stress,Displacements)

    Physics
     {
      stiffness isotropic
        {
          young = 149
          poisson = 0.349
        }

     }


On the upper surface we apply a normal force toward the base of the nanowire. The units are in GPa. ::

  Contact Upper 
  {
   type = surface_force
   force = (0,0,-0.00625)
  }


Finally, the boundary condition *clamp* frozes the column base. ::

  Contact Lower 
  {
   type = clamp
  }


As the simulation is at the equilibrium, we only need to solve the Poisson equation. This can be specified by :: 

  coupling = poisson  

The lower contact is ohmic. The output potential can be taken at the top of the nanowire, where the force is applied, with respect to the electrical potential at the base. 
This can be done by setting a reference contact, i.e. :: 

  reference_contact = Lower

With ::

  default_boundary_condition = zero_displacement

we enforce zero total displacement along the simulation boundary. 
Thus, here are   the options for  Module driftdiffusion ::


  Module driftdiffusion
  {

   coupling = poisson
   reference_contact = Lower
   default_boundary_condition = zero_displacement
   plot = (ElPotential,ElField,eDensity,Ec,Ev,Polarization,IonizedDonors)


By ::

  strain_simulation = elasticity

we specify the strain simulation, which will provide the strain map. ::

  Physics
  {
    strain_simulation = elasticity   
  }

    Contact Lower
    {
     type = ohmic
     regions = Lower
     voltage = 0.0
    }
  

To get the output piezopotential we first compute elasticity and then driftdiffusion. ::

  Simulation
   {
    solve = (elasticity,driftdiffusion)
    resultpath = output
    verbose = 3
   
    output_format = vtk
   }  


The Figure :ref:`Piezopotential` shows the electrical potential across a 200 :math:`nm` length region at the tip of the NW.
The output potential is :math:`\phi_M=-42mV`. This value is a result of the competition between the piezoelectric field and the screening due to the free carriers. 





.. _fig.piezopotential: 

.. figure:: ../data/figure_piezopotential.png
   :align: center
   :scale: 100%
   
   Piezopotential


   
..   </marker>




