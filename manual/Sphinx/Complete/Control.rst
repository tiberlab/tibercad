###########################
Controlling the Simulation
###########################

..
  .. index:: double:Simulations;sweep

================
Tools for Simulation
================


.. _Variable:


Variables
---------

In  tiberCAD,  the parameters defined in a module (e.g. voltage bias, temperature, material parameters) may be associated to a user-defined *variable*, in  such  a  way  that they their numerical value  can  be  controlled during the  simulation.
The  assignement to  a  *variable* is done  in  this  way ::

  parameter = $my_var[0.0]

where  ``my_var`` is  a  user-defined *variable* and the value in  square brackets [0.0] is  the  default  value of  the  *variable*.

For example ::

  Contact anode  
   {
    
    type = ohmic
    voltage = $Vbias[0.0]
   
   }

 
here the value of the  voltage in  the  defined ``Contact`` is  given  by the value of the *variable* ``Vbias``

 
A  typical application  of a user-defined *variable* is  in  the  *Sweep* calculation.
 


Sweep
---------

Performs a linear sweep for a given variable.


Options and  parameters:

 ``solve`` :  string
    a list of simulations to be solved, eg. ``(strain, dd)``

 ``variable`` : string
    the sweep variable, eg. ``$Vd``

 ``start`` :  double
    the first value of the sweep

 ``stop`` : double
    the last value of the sweep

 ``steps`` : integer
    the number of steps the interval (start, stop) gets subdivided in

 ``values`` :  double
  instead of ``start``, ``stop`` and ``steps``, provide the sweep values explicitly

 ``plot_data`` : boolean
    set to true if you want to plot after each step the simulation results (default is false)

 ``file_mode`` : string
    controls the behaviour for writing the data file containing global data. Can be one of ``append``, ``overwrite`` (default) or ``no-overwrite``

 ``min_step`` :  double
   the minimum absolute step size

 ``max_step`` : double
   the maximum absolute step size

 ``initial_step`` : double
   the absolute initial step size

 ``min_relative_step`` : double
     minimum relative step size, default is 1e-3

 ``max_relative_step`` : double
    maximum relative step size, default is 1

 ``initial_relative_step`` : double
    initial relative step size, default is 1

The relative step sizes refer to each single sweep step.  If  ``max_relative_step`` is less than one, each sweep step will be subdivided in smaller steps. If a simulation fails, the step gets reduced by half until the simulation succeeds, or until the minimum relative or absolute step size is reached. In the latter case, the sweep is assumed to have failed.
As for any module, a name can be given to a *sweep* block. This is important when several sweeps are defined, and in particular when nested sweeps (a solve option of one sweep refers to another sweep) are used.

Example::

  Module sweep
  {
   solve = (dd, thermal)
   variable = $Vbias
   start = 0.0
   stop = 1.2
   steps = 12
   plot_data = true
  }

In  this example, at  each  step of the  sweep the two simulations ``dd`` and  ``thermal`` are  performed, in  this  order, while  the  variable is ``Vbias`` and assumes values between 0 and 1.2.

 
..
  .. index:: double:Solvers;selfconsistent
	

.. _Selfcons:


Selfconsistent
----------------------

This tool may  be  used to solve  models in an iterative way to obtain a *selfconsistent* solution, optionally using a relaxation approach.

Options and parameters:

 ``name`` :  string
    optional name of  the ``selfconsistent``  block

 ``max_iterations`` : integer
    the maximum number of iterations. When the maximum number of iterations is reached, the program issues a warning, exits from the  self-consistent loop and proceeds to  the  next operation.

 ``relative_tolerance`` : double 
    the relative convergence tolerance for the observed variable in terms of the l2-norm

 ``absolute_tolerance`` : double
   the absolute convergence tolerance for the observed variable in terms of the maximum-norm

 ``relaxation_factor`` : double
   an optional relaxation factor (default is 1) to be applied to the observed variable

 ``solve`` :  (list of)  string(s)
   the simulations to be solved

The observed variable on which convergence control and relaxation is done is the system variable of the last simulation specified in  ``solve``.


Example::

  Module  selfconsistent
   {
     name = sc_all 
     solve = (quantum_el,quantum_hl, dd)
     max_iterations = 25  # 5

     absolute_tolerance = 1e-3
     relative_tolerance = 1e-5
   }

In this  example, the ``selfconsistent`` solution named ``sc_all`` is  defined,  which solves the  simulations defined in the relevant  modules as *quantum_el,quantum_hl* and *dd*,  in this  order. 


Then, in **Simulation**, we  execute ``sc_all``, e.g. ::

  Simulation
  {
   verbose = 3

   solve = (strain, sweep, scall)

   resultpath = output
   output_format = grace

}



=========================
Solvers
=========================

Numerical Solvers
-------------------------

The **Solver**  block inside a module description contains the options for the numerical solver.
The solver type the **Solver**  block is describing (*linear, nonlinear, eigenvalue solver*) depends on the module.
The options and  parameters for the *nonlinear* and the *linear* solvers are presented in  the following.

..
  .. index:: double:Solvers;nonlinear


.. _Nonlinear_solver:

Nonlinear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^

The Solvers implemented are of  two types :

  :ref:`petsc` =  uses the PETSc nonlinear solver (SNES) (default)

  :ref:`linesearch` =  uses a linear linesearch implemented in TiberCAD

Nonlinear solvers  are based on iterative methods, solving in each iteration a linear system. 
The linear solver used for this solution can be controlled by providing a block with the keyword ``linear_solver`` 
containing the options for the linear solver (see :ref:`Linear_solver`).



.. _petsc:



Petsc solver
...............


 ``relative_tolerance``	: double
    convergence criterion based on relative residual l_2-norm

 ``absolute_tolerance``	: double
    convergence criterion based on the l_2-norm of the residual

 ``max_iterations``:  integer
    maximum number of iterations

 ``step_tolerance`` : double
    tolerance criterion based on the l_2-norm of the correction step

 ``max_step`` :  double
    maximum linesearch step (l_2-norm)

 ``divergence_tolerance`` : double
    divergence criterion

..  _linesearch:


Linesearch solver
................


 ``absolute_tolerance`` : double
    convergence criterion based on the l_2-norm of the residual

 ``step_tolerance`` : double
    tolerance criterion based on the l_infinity-norm of the correction step

 ``max_iterations`` : integer
    maximum number of iterations

..
  .. index:: double:Solvers;linear


.. _Linear_solver:


Linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^


Petsc solver
...............



 ``method`` : string 
    The Krylov subspace method to be used. Here  are  the  available options:
  
 ``bcgs`` :
     BiCGstab

 ``bcgsl`` :
      BiCGstab(L)

 ``gmres`` :
    Generalized Minimal Residual

 ``bicg`` :
     BiConjugate Gradient

 ``cg`` :
     Conjugate Gradient

 ``cgs`` :
     Conjugate Gradient Squared

 ``richardson`` :
      Richardson

 ``pconly`` :
     only apply preconditioner

 ``relative_tolerance`` :
     convergence criterion based on relative residual l_2-norm

 ``absolute_tolerance`` :
    convergence criterion based on the l_2-norm of the residual

 ``max_iterations`` :
    maximum number of iterations


Available types of Preconditioner :

 ``lu``	:
    LU
 ``ilu`` :
    incomplete LU
 ``jacobi`` :
    Jacobi
 ``cholesky`` :
    Cholesky
 ``none`` :
     no preconditioning
 ``composite``	:
     use a combination of Jacobi and iLU


