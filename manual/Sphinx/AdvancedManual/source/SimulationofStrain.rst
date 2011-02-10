.. _SimulationofStrain:

Simulation of Strain
====================

Theory
------

The theoretical model of strain simulation can be found in [Cit-A]_ and [Cit-B]_ . 
The code can compute elastic deformations in a heterostructure and can calculate the deformed shape
of the structure. The heterostructure can be either grown on a substrate or not. External
pressure may be applied to a structure as well.

Models section parameters
-------------------------

The **Models** section looks like follows::

  model macrostrain
    {
     options
       {
        simulation_name = strain_in_transistor
        physical_regions = (2, 3, 4)
       }
     BC_Regions
       {
        .................
       }
    }
  
There are three possible kinds of boundary conditions. The mandatory keyword :
* type = { substrate | pressure | extended_material }
specifies the boundary condition type.

Substrate boundary condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case the boundary condition region ( see :ref:`DefinitionRegions` ) is the boundary between the device
and the substrate. The substrate does not belong to the device. Therefore it is necessary
to define both the boundary region number and the substrate material.

**Role of substrate**: In general, *the "substrate" is a material that defines the lattice
matching conditions, and not necessarily a real solid body on which the device is grown*::

  BC_Region layer_of_Al_0.3_Ga_0.7_N
    {
     BC_reg_numb = 14
     type = substrate
     mat = AlGaN
     x = 0.3
     structure = wz
     y-growth-direction = (1,0,-1,0)
     z-growth-direction = (-1,2,-1,0)
     x-growth-direction = (0,0,0,1)
    }
  
External pressure boundary condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameter *pressure* specifies the value (in GPa) of the normal pressure applied to
the boundary region *BC_reg_numb*::

  BC_Region tip_upon_a_surface
    {
     BC_reg_numb = 12
     type = pressure
     pressure = 12.3
    }
  
**Sign of pressure**. *The value of the pressure has a positive sign if the external force acts
towards the surface, which in general has to be the boundary of a simulation environment*.

Extended device boundary condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If device that is grown on a substrate is very large we may want to simulate a part of it
only. In this case the simulation domain boundary is not a free surface any more. The
boundary conditions are as follows:

.. math::
   :nowrap:
   :label:
   
   \begin{equation}
   \int C_{ijkl} \frac{\partial u_k}{\partial x_k} n_i d\partial\Omega = 0
   \end{equation}
   

The syntax is as follows::

  BC_Region boundary1
    {
     BC_reg_numb = 12
     type = extended_material
    }
  
Solver parameters
-----------------

The choice of the necessary parameters to be put in the **Solver** section depends on the 
type of the strain boundary condition for the structure, namely, if it is considered as 
grown on a substrate or not.

Structure with a substrate
^^^^^^^^^^^^^^^^^^^^^^^^^^

The only mandatory parameter is *substrate*, to which a name of a substrate boundary
condition region has to be assigned, e.g (referring to the previous example)::

  strain_in_transistor
    {
     substrate = layer_of_Al_0.3_Ga_0.7_N
    }

Structure without a substrate (freestanding)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, the parameter substrate should not be present; instead, the following pa
rameters should be defined.

* The reference lattice material is defined by the coordinates of a point belonging to
this material, using the parameter reference_material_point.

* As follows from [Cit-A]_ , additional geometrical points have to be specified, according 
to the device dimensionality. The corresponding parameters are: *fixed point1,
fixed point2 and fixed point3*.

Since the elasticity energy is invariant with respect to translations and rotations of the
structure, then, for the sake of uniqueness of solutions of the equations, another set of
constraints is required. Hereafter we assume that a mesh is defined over the simulation
domain and the displacement field :math:`u(r)` is defined at the mesh nodes.

Let *D* be the dimensionality of the structure minus the number of directions along
which the structure is periodic. If *D* > 0, then a Dirichlet boundary condition is applied
at an arbitrarily chosen node i1:

.. math::
   :nowrap:
   :label:
   
   \begin{equation} 
   {\bf u} ({\bf r}_{i_1}) = {\bf 0},
   \end{equation}

in order to prevent the structure from undesirable translations. In the case of D > 1,
another node :math:`i_2` is chosen and a constraint

.. math::
   :nowrap:
   :label: constr2D
   
   \begin{equation}
   ({\bf u} ({\bf r}_{i_1}) - {\bf u} ({\bf r}_{i_2}))\cdot({\bf r}_{i_1} - {\bf r}_{i_2}) = 0
   \end{equation}

is applied, in order to keep the direction between the nodes :math:`i_1, i_2` unchanged. If D = 3,
another node :math:`i_3` is chosen and additional constraint is set:

.. math::
   :nowrap:
   :label:
   
   \begin{equation}
   \left( {\bf r}_{i_1}- {\bf r}_{i_2}, {\bf r}_{i_1} - {\bf r}_{i_3}, {\bf u}_{i_3} \right) =0,
   \end{equation}

so the node i3 has to belong to the :math:`(\bf r_{i_1}, \bf r_{i_2}, \bf r_{i_3})` plane.

Example for a 2D simulation::

  strain_in_transistor
    {
     reference_material_point = (0, 100, 0)
     fixed_point1 = (0, 0, 0)
     fixed_point2 = (10, 0, 0)
    }

Additional parameters
^^^^^^^^^^^^^^^^^^^^^

**Numercal solver parameters**

*tolerance*: relative tolerance of the iterative solver, e.g.::

  tolerance = 1e-10

The default value is 1e-10.

*ksp type* = type of solver: **gmres, bcgsl, bcgs, cg, richardson**.
  Default is:

  * **gmres** for 1D
  * **bcgsl** for 2D and 3D
  
*pc* = type of pre-conditioner : ilu, composite, jacobi, lu, cholesky, eisenstat**.

Default is :

  * **ilu** for 1D
  * **jacobi** for 2D and 3D

*max_iterations* = max number of iterations (default = 1000)
*monitor, xmonitor* : if **true**, textual or graphical monitor of convergence process is
enabled (default = **false**)

**Periodic boundary conditions**

It is possible to specify periodic boundary conditions along the coordinate axes. The
relative parameters are:

    *periodicity_x* = { true | false }
    *periodicity_y* = { true | false }
    *periodicity_z* = { true | false }

    The default value is **false**.
  
**Mesh refinement**

For the details about mesh refinement refer to the Libmesh library documentation. The
parameter *refinement_steps* defines the number of the refinement steps to be done (with
default value equal to zero). The parameter *uniform_refinement* = { true | false } is
used to choose between uniform and adaptive refinement. The default value is **false**, i.e.
adaptive refinement. Example::

  refinement_steps = 4
  refine_fraction = 0.25
  coarsen_fraction = 0
  max_refinement_level = 10
  
**Deformed shape calculation**

The displacement field and lattice matching parameters that are found from master
equations can be used in order to define a new shape of the heterostructure. This new
shape is the first approximation to the equilibrium one. The next approximations are 
obtained iteratively by the following steps: at the n-th iteration the master equations
are solved using the lattice matching deformation :math:`\varepsilon^n_{ij}` which is defined as

.. math::
   :nowrap:
   :label: eps1
   
   \begin{equation}
   \varepsilon_{ij}^n = \frac{1}{2}\left(\frac{\partial u_i^{n-1}}{\partial x_j} +
                               \frac{\partial u_j^{n-1}}{\partial x_i}\right)
   + \epsilon_{ij}^{n-1},
   \end{equation}

where the displacement field :math:`u^{n-1}` has been taken from the iteration n - 1. Then the
new shape is defined by using the displacements from the last step solution, and the
iterative process is repeated until the displacement field vanishes and additional lattice
parameters stabilise. The iterative cycle usually converges after 3 - 4 iterations.

The only parameter that controls shape calculation is number_shape_steps. The value
defines number of iterations. The default value is zero, that means no shape deformation
calculation.

Physics section parameters
------------------------------

There is a possibility to consider converse piezoelectric effect. For this it is necessary to
specify a name of another simulation that can provide electric field. The parameter is
*poisson_equation*. Example::

  macrostrain
    {
     poisson_equation = DriftDiffusion
    }

Interaction with other simulations. *In order to take into account the converse piezo
effect, the poisson equation has to recalculate the necessary parameters after the strain
simulation. To do so, the following parameters has to be set in the Physics section of the
drift-diffusion equation* (for detailes see :ref:`PhysicSection` ) ::

  driftdiffusion
    {
     model = strained
     strain_simulation = str
     recompute_band_parameters = true
    }
  
Output
----------

The output variables are:

* strain strain tensor (6 components) in calculation system
* stress stress tensor (6 components) in calculation system
* polarization piezo polarization vector (3 components) in calculation system


  +---------------+--------------+-------------+
  | Plot keyword  | label        | Units       |
  +---------------+--------------+-------------+
  | strain        | eps          |             |
  +---------------+--------------+-------------+
  | stress stress | GPa          |             |
  +---------------+--------------+-------------+
  | polarization  | polarization | |Table5.1A| |
  +---------------+--------------+-------------+ 

Table 5.1: Elemental vector quantities

By using **StrainVariables** as a plot keyword it's possible to include all quantities
of the strain simulation.


.. |Table5.1A| replace:: :math:`C/cm^{-2}`

.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes
