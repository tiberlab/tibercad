.. _EFA:



Module efaschroedinger
=================================


The EFA calculations are performed by the **Module** efaschroedinger A typical ex-
ample is the following::

  Module efaschroedinger
    {
     particle = el
     poisson_model_name = dd
     strain_model_name = strain # macrostrain
     name = quantum_el
     regions = quantum
     plot = (EigenFunctions, EigenEnergy, EnergyLevels)
     Solver
       {
        number_of_eigenstates = 10 # 30 
       }
     Physics
       {
        model = conduction_band
       }
    }

Module options
-----------------------

The following options infiuence the behaviour of the Module efaschroedinger:

  ``particle = string`` defines for which particle (electron or hole) Schroedinger equation is
solved Possible values are el and hl. A different Module efaschroedinger has to be
defined for each particle to be solved.

  ``poisson_model_name = string`` defines the name of the simulation (Module driftdiffu-
sion) that can provide electric potential

  ``strain_model_name = string`` defines the name of the simulation (Module macrostrain)
that can provide elastic strain

  ``regions = string`` defines the regions associated to this EFA simulation

Solver section
--------------------------

The Solver section of the Module efaschroedinger contains the following options:

  ``number_of_eigenstates = integer`` defines the number of eigenvalues and eigenfunctions to be found.

  ``Dirichlet_bc_everywhere = boolean`` : if true (default value), Dirichlet boundary
conditions are imposed over all the boundaries of the simulation region

  ``solver = string`` : defines the solver for the eigenvalue problem, possible values are:
arnoldi, lapack, krylovshur. The default value is **krylovshur** . 
In the case of the lapack solver all the eigenvalues are computed. In the case of arnoldi or krylovshur
solver it is necessary to specify which and how many eigenvalues have to be computed. 
The idea is that the iterative solver calculates several eigenvalues that are
close to a specific number, refered to as the *spectral_shift*

  ``max_iteration_number = integer`` : maximum number of iteration, used as a stop condition

  ``eigen_solver_tolerance = double`` :numerical eigensolver tolerance used as a convergence criteria

  ``spectral_shift = double`` :the closest eigenvalues to this value (eV) are found. If not
defined, then it will be calculated internally from the band edges.

  ``ksp_type = string`` : Krylov subspace method type; bcgsl, gmres, cg

  ``pc_type = string`` : preconditioner type: cholesky, jacobi, ilu , composite.

Physics section
---------------------

  ``model = string`` : possible values are : conduction band , for single conduction band
model ( :math:`\Gamma` point) ;  kp for :math:`{\bf k \cdot p}` model

  ``kp_model = string`` : possible values are : 6×6 , 8×8.

Output
----------------------

The available output variables, to be specified in the plot option, are the following:

* EigenEnergy Eigen energy in eV

* EigenFunctions :math:`|\psi({\bf r})|^2` function of the eigenstate

* Occupation probability to find the state occupied. It is calculated assuming Fermi
distribution and mean electrochemical potential and temperature:

..  math::
    :nowrap:
	:label:

    \begin{gather}
    \bar \mu = \langle \psi|\mu({\bf r})|\psi \rangle \\
    \bar T = \langle \psi|T({\bf r})|\psi \rangle 
    \end{gather}

* EnergyLevels graphical output used for showing the energy level over the band
diagram.

Module opticskp
=========================

By defining the Module **opticskp** , 
calculation of optical properties is enabled.

::

  Module opticskp
    {
     name = optics
     regions = quantum
     plot = (optical_spectrum_k_0 )
     initial_state_model = quantum_el
     final_state_model = quantum_hl
     initial_eigenstates = (0, 9)
     final_eigenstates = (0, 15)
     polarization = (0, 0, 1)
     Emin = 2.8
     Emax = 3.6
     dE = 0.001
    }

Here, *initial_state_model* and *final_state_model* are, respectively, the quantum simulations 
( **efaschroedinger** module) associated respectively to the initial state of optical
transition (e.g. electron), and to the final state of optical transition (e.g. hole). 
*initial_eigenstates* and *final_eigenstates* refer to the range of eigenstates to be taken in
account for optical calculations.

By specifying a range of energy values in this way::

  Emin = 3.0
  Emax = 5.0
  dE = 0.001

the emission optical spectrum for **k=0** is calculated.

Output
---------------------------

The output variables for optics calculations are:

* **optical_spectrum_k_0** : optical emission spectrum for *k=0*.

Module opticalspectrum
================================

By defining the Module **opticalspectrum** , optical matrix elements are used to calculate
the associated (emission) spectrum with a k-space integration.

::

  Module opticalspectrum
    {
     k_space_dimension = 2
     k-space_basis = true
     k1 = (0, 0, 0.1)
     k2 = (0, 0.1, 0)
     refine_fraction = 0.30
     relative_accuracy = 0.01
     refine_k_space = true
     number_of_nodes = (2, 2)
     wedge = quarter
     plot = (optical_spectrum)
     optical_matr_elem_model = opticskp
     polarization = (0, 0, 1)
     Emin = 3.0
     Emax = 5.0
     dE = 0.001
    }

The parameters are the following:

  ``k_space_dimension`` = **1** for 2D simulations, **2** for 1D simulations. k-space basis is
 **true** if the k-space is defined by means of k-vectors; if **false** , vectors are expressed in
real space

  If ``refine_k_space`` = **true** , that is adaptive k-mesh refinement is enabled, all the el-
ements whose error is greater than the value (1-refine fraction)* (maximum error) are
going to be refined. In this case, ’Error’ is just the integrated quantity. The refinement
will end when the *relative_accuracy* is obtained.

  ``number_of_nodes`` = numb. of elements in k mesh, along each direction

  ``wedge`` = half | quarter, to reduce calculation time, by exploiting symmetry.

  ``optical_matr_elem_model`` = name of the *opticskp* model associated

  ``polarization`` = light polarization (vector)

  ``Emin, Emax, dE`` : energy range and step of spectrum calculation.

Output
---------------------

The output variables for optics calculations are:

* **optical_spectrum** : k-space integrated optical emission spectrum.Chapter 5

Module quantumdispersion
=============================

With the Module quantumdispersion it is possible to calculate the dependence of quan-
tum eigenstates on **k** -vector. Such dependence gives the *quantum state dispersion* . The
simulation name is **quantumdispersion** .

::

  Module quantumdispersion
    {
     simulation_name = dispersion1D_el
     regions = all
     quantum_simulation = quantum_el
     min_eigenvalue_number = 0
     max_eigenvalue_number = 5
     wedge = half
     k_space_dimension = 1
     k1 = (0, 0.1, 0)
     number_of_nodes = (10)
     output_format = grace
     plot = k-space_dispersion
    }

The dispersion of quantum states is calculated at k-points that are nodes of the mesh
in k-space.

The main parameters are:

* **quantum simulation** : name of the efaschroedinger simulation.

* **min eigenvalue number** , **max eigenvalue number** : the dispersion is calculated 
for the states number *i* , where 

  max_eigenvalue_number :math:`\ge i \ge`  min_eigenvalue_number

The rest of the parameters (wedge, k space dimension, etc...) define the k-space.

Output
--------------------

The output variable name is **k-space_dispersion** . The output format for the dispersion
can be controlled independently of the general specification in the ``Simulation`` section
by redefining the ``output_format`` keyword.

Module quantumdensity
==============================

In Module quantumdensity you can define the calculation of particle (electron,hole) den-
sity, based on the result of a previous calculation of the system eigenstates (e.g. with
efaschroedinger module). 

:: 

  Module quantumdensity
    {
     name = dens_el
     regions = quantum
     plot = quantum_density
     k-space_dimension = 0
     k-space_basis = true
     k1 = (0, 0, 0.1)
     #number_of_nodes = (4)
     #wedge = half
     quantum_simulation = quantum_el
     degeneracy = 2
     refine_fraction = 0.20
     relative_accuracy = 0.01
     refine_k_space = true
     output_density_in_k_space = true
     uniform_refinement = false
     mesh_order = FIRST
     first_state = 3
     initial_eigenstates_number = 10
     analytic = true
    }

The available options are:

* **quantum_simulation** : name of the efaschroedinger simulation.

* degeneracy: degeneracy of the quantum state

* **initial_eigenstates_number** : initial number of eigenstates for the Schroedinger
equation

* **analytic** = { true | false } If true then the density is calculated analytically,
otherwise numerically.

Output
====================

The output parameter is **quantum_density**.





.. |more| image:: ../data/more.*
    :scale: 50%

.. |warn| image:: ../data/warn.*
    :scale: 50%

.. |idea| image:: ../data/idea.*
    :scale: 50%


.. rubric:: Footnotes

