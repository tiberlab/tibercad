.. _EnvelopeFunction:

Envelope Function Approximation
===============================

The envelope function approximation (EFA) simulation tool of TIBERCAD is developed
in order to solve a single-particle Schroedinger equation for electrons and holes in a semi
conductor crystal. This problem is an eigenvalue problem that is treated as a generalized
complex eigenvalue problem

.. math::
   :nowrap:
   :label:
   
   \[
   H\psi = ES\psi 
   \]

where H and S are the Hamiltonian and S-matrix, respectively.

Models section parameters
-------------------------

The **Models** section looks like follows::

  model efaschroedinger
    {
     options
       {
        simulation_name = quantum_well1
        physical_regions = (1,2)
       }
    }
  
The default boundary conditions of the simulation domain are open (that is zero flux
for single-band calculation). It is possible to specify Dirichlet boundary conditions::

  BC_Region infinite_barrier1
    {
     BC_reg_numb = 12
     type = Dirichlet
    }

There is a way to impose automatically the Dirichlet boundary conditions over all
the boundary of the simulation region. This is done by the parameter:

*Dirichlet_bc_everywhere* = {true | false} in **Solver** section. The default value for EFA
problem is true.

Solver parameters
-----------------

There are two groups of parameters. The first group is related to the general eigensolver
problem, the second one is related to the Schroedinger equation.

Eigenvalue problem parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters are common for all eigenvalue problems. Their default values may be
different for different eigenvalue problems, for example for the Schroedinger equation and
for the electromagnetic eigenvalue problem.

The eigenvalue problem can be solved by the solvers that are implemented into the
SLEPc library. The relative parameter is

*solver* = { arnoldi | lapack | krylovshur }. The default value is **krylovshur** .

In the case of the lapack solver all the eigenvalues are computed. In the case of arnoldi
or krylovshur solver it is necessary to specify which and how many eigenvalues have to
be computed. 
The idea is that the iterative solver calculates several eigenvalues that a
close to a specific number, referred to as the **spectral shift** . The relative parameters are:

| 

+----------------------------------+---------------------------------------------------------------+
| attribute                        | description                                                   |
+----------------------------------+---------------------------------------------------------------+
| max_iteration_number             | maximum number of iteration, used as a stop condition         |
+----------------------------------+---------------------------------------------------------------+
| eigen_solver_tolerance           | numerical eigensolver tolerance used as a convergence criteria|
+----------------------------------+---------------------------------------------------------------+
| spectral_shift                   | the closest eigenvalues to this value (eV) are found.         |
+----------------------------------+---------------------------------------------------------------+
| spectrum_inversion               | tolerance tolerance used for linear solver                    |
+----------------------------------+---------------------------------------------------------------+

Table 7.1: Iterative eigensolver parameters

| 

If the spectral shift is not defined then it will be calculated internally from the band
edges.

For the iterative solvers the important parameters, that may significantly change the
performance, are the *Krylov subspace method type* and the *preconditioner* type.
The Krylov method is defined as follows:
ksp type = { bcgsl | gmres | bcgs | cg | richardson | preonly }

The **preconditioner** type is defined as follows:

pc_type = { cholesky | jacobi | ilu | composite }

Other options:

* *x-periodicity* = { true | **false** }

* *y-periodicity* = { true | **false** }

* *z-periodicity* = { true | **false** }

* *number of eigenstates* = number of eigenstates to be computed

Schrodinger equation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* *particle* = { el | hl }; specifies a particle type, according to each the eigenvalues are sorted

* *strain_model name* = name of the simulation that can provide elastic strain

* *poisson_model name* = name of the simulation that can provide electric and electrochemical potential

* *heat_model* = name of the simulation that can provide

* *relative_density tolerance* = relative tolerance for charge density calculation

* *initial_eigenstates number* = initial eigenstates number for charge density calculation

* *convergent_density* = { true | **false** } if true than number of eigenstates will be
increased in oder to reach convergent density.

* eigen_number_increase_factor = factor to increase eigenvalues number for the next
charge density calculation

Physical Models parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

* *particle* = {el | hl }

* *model* = { conduction band | kp } single conduction band ( point) or **:math:`k \cdot p`** .

* *kp_model* = { 6x6 | 8x8 }

Here, the *particle* name is the name of a particle type (*electron* or *hole*).

*model* = **kp | conduction band** : **:math:`k \cdot p`** or single conduction band model.
If **:math:`k \cdot p`** model is applied, specify::

  *kp_model* = **6x6 | 8x8**.
  
Output
------

.. index:: single:Eigenstatus

* **EigenEnergy** Eigen energy in eV

* **EigenFunctions** :math:`|\phi(r)^2|` function of the eigenstate

* **Occupation** probability to find the state occupied. It is calculated assuming Fermi
distribution and mean electrochemical potential and temperature:

.. math::
   :nowrap:
   :label:
   
   \begin{align*}
   \mu =  \left  \langle \psi|\mu (\mathbf{r})|\psi  \right \rangle \\
   T = \left  \langle \psi| T (\mathbf{r})|\psi  \right \rangle \\
   \end{align*}
   
* **EnergyLevels** graphical output used for showing the energy level over the band
diagram





.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes