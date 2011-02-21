.. _DriftDiffusion:

Solvers
=========================

Numerical Solvers
-------------------------

the “Solver” block inside a module description contains the options for the numerical solver.
The solver type the “Solver” block is describing (linear, nonlinear, eigenvalue solver) depends on the module.
For nonlinear and linear solvers, the following options exist:

.. index:: double:Solvers;nonlinear

Nonlinear solvers:
^^^^^^^^^^^^^^^^^^^^^^^^^

type :

  :ref:`petsc` =  uses the PETSc nonlinear solver (SNES) (default)

  :ref:`linesearch` =  uses a linear linesearch implemented in TiberCAD

Nonlinear solvers  are based on iterative methods, solving in each iteration a linear system. 
The linear solver used for this can be controlled by providing a block with keyword “linear_solver” 
containing the options for the linear solver (see linear solvers).

.. _petsc:

**petsc**

* ``relative_tolerance``	=  convergence criterion based on relative residual l_2-norm

* ``absolute_tolerance``	=  convergence criterion based on the l_2-norm of the residual

* ``max_iterations``	=  maximum number of iterations

* ``step_tolerance`` =  tolerance criterion based on the l_2-norm of the correction step

* ``max_step`` =  maximum linesearch step (l_2-norm)

* ``divergence_tolerance	=  divergence criterion

..  _linesearch:

**linesearch**

* ``absolute_tolerance`` =  convergence criterion based on the l_2-norm of the residual

* ``step_tolerance`` =  tolerance criterion based on the l\infinity-norm of the correction step

* ``max_iterations`` =  maximum number of iterations

.. index:: double:Solvers;linear

Linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^

type = petsc

**Petsc**

  ``method``, The Krylov subspace method to be used:
  
  ``bcgs`` =  BiCGSTAB
  ``bcgsl`` =
  ``gmres``	=  Generalized Minimal Residual
  ``bicg`` =  BiConjugate Gradient
  ``cg`` =  Conjugate Gradient
  ``cgs`` =  Conjugate Gradient Squared
  ``richardson`` =  Richardson
  ``pconly`` = only apply preconditioner

* ``relative_tolerance`` =  convergence criterion based on relative residual l_2-norm

* ``absolute_tolerance`` =  convergence criterion based on the l_2-norm of the residual

* ``max_iterations`` =  maximum number of iterations

Preconditioner :

  ``lu``	=  LU
  ``ilu``	=  incomplete LU
  ``jacobi`` =  Jacobi
  ``cholesky`` =  Cholesky
  ``none``	=  no preconditioning
  ``composite``	=  use a combination of Jacobi and iLU

.. index:: double:Simulations;sweep

Simulations
---------------------

**Sweep**

Performs a linear sweep for a given variable.

Options:

  ``solve`` =  a list of simulations to be solved, eg. (strain, dd)

  ``variable`` =  the sweep variable, eg. $Vd

  ``start`` =  the first value of the sweep

  ``stop`` =  the last value of the sweep

  ``steps`` =  the number of steps the interval (start, stop) gets subdivided in

  ``values`` =  instead of start, stop and steps, provide the sweep values explicitly

  ``plot_data`` =  set to true if you want to plot after each step the simulation results (default is false)

  ``min_step`` =  the minimum absolute step size

  ``max_step`` = the maximum absolute step size

  ``initial_step`` = the absolute initial step size

  ``min_relative_step`` =  minimum relative step size, default is 1e-3

  ``max_relative_step`` =  maximum relative step size, default is 1

  ``initial_relative_step`` =  initial relative step size, default is 1

The relative step sizes refer to each single sweep step.  If max_relative_step is less than one, each sweep step will be subdivided in smaller steps. If a simulation fails, the step gets reduced by half until the simulation succeeds, or until the minimum relative or absolute step size is reached. In the latter case, the sweep is assumed to have failed.
As for any module, a name can be given to a sweep block. This is important when several sweeps are defined, and in particular when nested sweep (solve option of one sweep refers to another sweep) are used.

Example::

  Module sweep
    {
     solve = (dd, thermal)
     start = 0
     stop = 1
     steps = 10
     plot_data = true
    }

.. index:: double:Solvers;selfconsistent
	
Selfconsistent
----------------------

Solves models in an iterative way to obtain a selfconsistent solution, optionally using a relaxation approach.
options:

* ``max_iteration`` =  the maximum number of iterations. Currently, when the maximum number of iterations is reached, the program only issues a warning and proceeds.

* ``relative_tolerance`` =  the relative convergence tolerance for the observed variable in terms of the l2-norm

* ``absolute_tolerance`` =  the absolute convergence tolerance for the observed variable in terms of the maximum-norm

* ``relaxation_factor`` =  an optional relaxation factor (default is 1) to be applied to the observed variable

* ``solve`` =  the simulations to be solved

Currently, the observed variable on which convergence control and relaxation is done is the system variable of the last simulation specified in 'solve'.




