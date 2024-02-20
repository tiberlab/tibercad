###########################
Controlling the Simulation
###########################

..
  .. index:: double:Simulations;sweep

================
Tools for Simulation
================

.. _Parallel:


Parallel execution
---------

tiberCAD may be executed in parallel on a number *n* of processors by typing:

   
  ``tibercad`` -n   **<number_MPI_proc>**  *input_file_name*


where **<number_MPI_proc>** is the chosen number *n* of processors.

In the ``Device`` section it is possible to define a set of options for the parallel execution 
using a block with the keyword ``Parallel``.

For example:

::
   
  Parallel
  {
    mpi_processes_per_device = 2 # 4 
    
    mpi_processes_per_mesh   = 1
  }


The  following options may  be  defined in **Parallel**  block


 ``mpi_processes_per_device`` : integer, (*default = number_MPI_proc*)       
    defines the number of devices to be processed in parallel, given by *number_MPI_proc/mpi_processes_per_device* 

 ``mpi_processes_per_mesh`` : integer , (*default = mpi_processes_per_device*)
    defines the number of processes for FEM calculations; if 1, then FEM calculations are executed in serial; 
    max value is *mpi_processes_per_device*

By default, FEM solvers will be parallelized on the defined number of MPI processes, that is ``mpi_processes_per_device`` and  ``mpi_processes_per_mesh`` are set to the value of **number_MPI_proc**.

All the other solvers, such as atomistic solvers, will be  parallelized on the defined number of MPI processes.


When a simulation requires a calculation in **k**-space, such as a **k**-integration, this calculation is automatically parallelized on a number of processes given by 
*mpi_processes_per_device/mpi_processes_per_mesh*

Current limitations
^^^^^^^^^^^^^^^^^^^^^^^

It is not possible to perform a parallel simulation of a device where different portions of the mesh are assigned to  different Modules. 


.. warning::
             A Simulation can be parallelized only if **ALL** the Modules are applied to the same set of Physical Regions


It is not currently possible to run Module DSC in parallel.


Examples
^^^^^^^^^^^^^^^^^^^^^^^

In multiscale simulations including both FEM-based and atomistic Modules (e.g. ETB), it is convenient to run in parallel only the atomistic solver, much more time-consuming and computational heavy, while FEM Modules, such as Drift-Diffusion, can be  solved serially.
 
This can be done quite easily in this way:

::

	Device
	{
  	meshfile = InAs_qw.msh 

  	# if run in parallel, only TB part will be parallelized
  	Parallel
  	{
    	    mpi_processes_per_mesh = 1
 	}
 	...........................

Here, when executed in parallel, only the atomistic solvers (**VFF**, **ETB**) will be parallelized.


When running in parallel with MPI, the default is to have one single device and run solvers in parallel.
It may be the case that  we want to run simulations of a certain number of devices in parallel, with or without parallel solvers.
This can happen for example when one wants to perform a set of simulations of a device including a random alloy material in the active region, so that to  generate a statistical set of solutions.
This is possible by defining appropriately the keyword *mpi_processes_per_device*.

For example with 

**number_MPI_proc** = 4 

let's consider a Tight-binding calculation on an atomistic random alloy structure: 

::

  	Parallel
  	{
    	    	mpi_processes_per_device = 2 
		# let's do FEM calculations in serial
    		mpi_processes_per_mesh   = 1

 	}


        ...........................


 	Atomistic tb
  	{
    	reference_region = qbarrierl
    	regions = atomistic
    	passivation = yes
    	print = (xyz)

    	random_alloy = true

    	# we put an implicit seed here, for repeatability (and reload)
   	 random_generator_seed = MPI_DEV_KEY
	# same random seed for each device process group

    	supercell_size_z = 20 
    	supercell_size_y = 20  

    	# write out statistics of alloy
    	alloy_statistics = true
    	alloy_statistics {}

  	}


Two devices  will be simulated (with a different random alloy), each one with a parallelization on 2 processes.
This is set by ::

   mpi_processes_per_device = 2



Note that in this case *mpi_processes_per_mesh* may be only 1 (serial) or 2 (FEM paralleized on 2 processes, which is the default in this case)
If *mpi_processes_per_device = 1* , then 4 different devices are simulated, with no parallelization.
In this case *mpi_processes_per_mesh* must be 1.

If *mpi_processes_per_device = 4*, then only one device is simulated, with a parallelization on 4 processes


Note that with *random_generator_seed* we fix the same random seed for each device process group, so that parallelized jobs keep the same random alloy  structure.
 


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

Performs a linear or logarithmic sweep for a given variable.



Options and  parameters:

 ``type`` :  string
    type of sweep, ``linear`` (default) or ``log``

 ``solve`` :  string
    a list of simulations to be solved, eg. ``(strain, dd)``

 ``variable`` : string
    the sweep variable, eg. ``$Vd``

 ``start`` :  double
    the first value of the sweep

 ``stop`` : double
    the last value of the sweep

 ``steps`` : integer
    the number of steps the interval (start, stop) gets subdivided in.
    For logarithmic sweeps, default number of steps is calculated assuming
    a factor of 10 between successive values. 

 ``values`` :  double
  instead of ``start``, ``stop`` and ``steps``, provide the sweep values explicitly. This accepts also ranges in Matlab syntax, e.g. ``0:0.25:4``

 ``plot_data`` : string
    set to ``each`` if you want to plot after each step the simulation results (default is ``none``), or to ``last`` to write only results of last step. Alternatively, a vector of values can be given, indicating the specific values where to plot data. Values not present in ``values`` or the variable range will be ignored.

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
   plot_data = (0, 0.3, 0.7, last)
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

 ``convergence_check`` : string [optional]
   the name of the simulation to be used for checking convergence

The observed variable on which convergence control and relaxation is done is the system variable of the last simulation specified in  ``solve``, if
not specified explicitly by the ``convergence_check`` option.


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
    convergence criterion based on relative residual l_2-norm, default 1e-6

 ``absolute_tolerance``	: double
    convergence criterion based on the l_2-norm of the residual, default 1e-50

 ``max_iterations``:  integer
    maximum number of iterations, default 500

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
    convergence criterion based on the l_2-norm of the residual, default 1e-50

 ``relative_tolerance``	: double
    convergence criterion based on relative residual l_2-norm, default 1e-9 

 ``step_tolerance`` : double
    tolerance criterion based on the l_infinity-norm of the correction step , default 1e-3

 ``max_iterations`` : integer
    maximum number of iterations , default 25

..
  .. index:: double:Solvers;linear


.. _Linear_solver:


Linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^


Petsc solver
...............

The Krylov subspace method to be used can be selected with the keyword **method**. 
For example ::
  
  method = cgs   


Here  are  the  available values for *method*; default is **bcgs** :
  
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


Other options :

 ``pconly`` :
     only apply preconditioner

 ``relative_tolerance`` :
     convergence criterion based on relative residual l_2-norm, default 1e-6 

 ``absolute_tolerance`` :
    convergence criterion based on the l_2-norm of the residual, default 1e-50 

 ``max_iterations`` :
    maximum number of iterations, default 1500


Preconditioner can be selected with the keyword **preconditioner**. 
For example ::
  
  preconditioner = ilu   

Available types of Preconditioner are :

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


