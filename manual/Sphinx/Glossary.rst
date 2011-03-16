.. _Glossary:

Glossary
========

.. glossary::
   :sorted:
   
   meshfile
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   meshfile = bulk.msh
      | 
      
   Doping
      | Level: sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |       Doping
      |       {
      |         Nd = 1e16
      |         type = donor
      |       }
      |
      
   Nd
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Nd = 1e16
      | 
      
   type
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   type = donor
      |   type = acceptor
      |   type = field_dependent
      |   type = heat_reservoir
      |   type = surface_resistance
      |   type = surface_force
      | 
      
   plot
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   plot = (Ec, Ev, Eg, ElField, ElPotential, eQFermi, hQFermi, eDensity, hDensity,
      |           NetRecombination, eCurrentDensity, hCurrentDensity, CurrentDensity,
      |           ContactCurrent, eMobility, hMobility, IonizedDonors, IonizedAcceptors)
      |   plot = all
      |
      
   solve
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   solve = driftdiffusion
      |   solve = sweep
      | 
      
   variable
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   variable = $Vb
      | 
      
   start
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   start = 0.0
      | 
      
   stop
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   stop = 1
      | 
      
   steps
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   steps = 10
      | 
      
   plot_data
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   plot_data = true
      | 
      
   verbose
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   verbose = 2
      | 
      
   resultpath
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   resultpath = output
      | 
      
   output_format
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   output_format = grace
      | 
      
   material
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   material = Si
      | 
      
   name
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   name = dd
      | 
      
   statistics
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   statistics = boltzmann
      |   statistics = FermiDirac
      | 
      
   low_field_model
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   low_field_model = doping_dependent 
      | 
      
   recombination
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   recombination (srh, auger) {}
      | 
      
   structure
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   structure = wz
      | 
      
   x-growth-direction
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   x-growth-direction = (0,0,0,1)
      | 
      
   y-growth-direction
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   y-growth-direction = (1,0,-1,0)
      | 
      
   z-growth-direction
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   z-growth-direction = (-1,2,-1,0)
      | 
      
   density
      | Level: paragraph
      | 
      | Description: doping concentration [cm-3]
      | 
      | 
      | Example: 
      |   density = 1e17
      | 
      
   level
      | Level: paragraph
      | 
      | Description: energy level of the dopant [eV]
      | 
      | 
      | Example: 
      |   level = 0.025
      | 
      
   strain_simulation
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   strain_simulation = macrostrain
      |   strain_simulation = elasticity
      | 
      
   nonlinear_solver
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   nonlinear_solver = tiber
      | 
      
   nonlin_max_it
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   nonlin_max_it = 15
      | 
      
   nonlin_rel_tol
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   nonlin_rel_tol = 1e-12
      | 
      
   nonlin_abs_tol
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   nonlin_abs_tol = 1e-10
      | 
      
   nonlin_step_tol
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   nonlin_step_tol = 1e-5
      | 
      
   polarization
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   polarization (piezo, pyro) {}
      | 
      
   max_iterations
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   max_iterations = 5
      | 
      
   abs_tolerance
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   abs_tolerance = 1e-3
      | 
      
   rel_tolerance
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   rel_tolerance = 1e-6
      | 
      
   mesh_units
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   mesh_units = 1e-6
      | 
      
   x
      | Level: paragraph
      | 
      | Description: the elementum percentage in the alloy
      | 
      | 
      | Example: 
      |           material = AlGaN
      |           x = 0.14
      | 
      
   H
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |           H = 1e10
      | 
      
   mesh_regions
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   mesh_regions = (barrier_dop_s, barrier_dop_d)
      | 
      
   regions
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   regions = all #default to not specify
      |   regions = -passivation #apply to all regions except passivation
      | 
      
   coupling
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   coupling = electrons
      | 
      
   save_state
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   save_state = true
      | 
      
   relative_tolerance
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   relative_tolerance = 1e-15
      | 
      
   step_tolerance
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   step_tolerance = 5e-3
      | 
      
   ksp_type
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   ksp_type = pconly #preconditioner only
      | 
      
   preconditioner
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   preconditioner = lu
      | 
      
   region
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   region = GaN/SiN
      |   region = buffer
      | 
      
   reference
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   reference = cb
      | 
      
   dimension
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   dimension = 3
      | 
      
   reference_material
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   reference_material = AlGaN
      | 
      
   temperature
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   temperature = 300
      | 
      
   r_surf
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   r_surf = 0. 05
      | 
      
   dd_simulation
      | Level: paragraph
      | 
      | Description: The name of the drift-diffusion simulation is given by dd_simulation.
      | 
      | 
      | Example: 
      |   dd_simulation = dd
      | 
      
   particle
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   particle = el
      | 
      
   poisson_model_name
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   poisson_model_name = dd
      | 
      
   strain_model_name
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   strain_model_name = strain # macrostrain
      | 
      
   number_of_eigenstates
      | Level: paragraph
      | 
      | Description: defines the number of eigenvalues and eigenfunctions to be found.
      | 
      | 
      | Example: 
      |   number_of_eigenstates = 10 # 30
      | 
      
   model
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   model = conduction_band
      | 
      
   simulation_name
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   simulation_name = dispersion1D_el
      | 
      
   quantum_simulation
      | Level: paragraph
      | 
      | Description: name of the efaschroedinger simulation.
      | 
      | 
      | Example: 
      |   quantum_simulation = quantum_el
      | 
      
   min_eigenvalue_number
      | Level: paragraph
      | 
      | Description: the dispersion  calculated for the states number *i*
      | 
      | 
      | Example: 
      |   min_eigenvalue_number = 0
      | 
      
   max_eigenvalue_number
      | Level: paragraph
      | 
      | Description: the dispersion  calculated for the states number *i*
      | 
      | 
      | Example: 
      |   max_eigenvalue_number = 5
      | 
      
   wedge
      | Level: paragraph
      | 
      | Description: to reduce calculation time, by exploiting symmetry.
      | 
      | 
      | Example: 
      |   wedge = half
      |   wedge = quarter
      | 
      
   k1
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   k1 = (0, 0, 0.1)
      | 
      
   k2
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   k2 = (0, 0.1, 0)
      | 
      
   number_of_nodes
      | Level: paragraph
      | 
      | Description: number of elements in k mesh, along each direction.
      | 
      | 
      | Example: 
      |   number_of_nodes = (10)
      |   number_of_nodes = (2, 2)
      | 
      
   k-space_dimension
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   k-space_dimension = 0
      | 
   k-space_basis
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   k-space_basis = true
      | 
      
   degeneracy
      | Level: paragraph
      | 
      | Description: degeneracy of the quantum state
      | 
      | 
      | Example: 
      |   degeneracy = 2
      | 
      
   refine_fraction
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   refine_fraction = 0.20
      | 
      
   relative_accuracy
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   relative_accuracy = 0.01
      | 
      
   refine_k_space
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   refine_k_space = true
      | 
      
   output_density_in_k_space
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   output_density_in_k_space = true
      | 
      
   uniform_refinement
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   uniform_refinement = false
      | 
      
   mesh_order
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   mesh_order = FIRST
      | 
      
   first_state
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   first_state = 3
      | 
      
   initial_eigenstates_number
      | Level: paragraph
      | 
      | Description: initial number of eigenstates for the Schroedinger equation
      | 
      | 
      | Example: 
      |   initial_eigenstates_number = 10
      | 
      
   analytic
      | Level: paragraph
      | 
      | Description: { true | false } If true then the density is calculated analytically, otherwise numerically.
      | 
      | 
      | Example: 
      |   analytic = true
      | 
      
   initial_state_model
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   initial_state_model = quantum_el
      | 
      
   final_state_model
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   final_state_model = quantum_hl
      | 
      
   initial_eigenstates
      | Level: paragraph
      | 
      | Description: range of eigenstates to be taken in account for optical calculations.
      | 
      | 
      | Example: 
      |   initial_eigenstates = (0, 9)
      | 
      
   final_eigenstates
      | Level: paragraph
      | 
      | Description: range of eigenstates to be taken in account for optical calculations.
      | 
      | 
      | Example: 
      |   final_eigenstates = (0, 15)
      | 
      
   polarization
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   polarization = (0, 0, 1)
      | 
      
   Emin
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Emin = 2.8
      | 
      
   Emax
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Emax = 3.6
      | 
      
   dE
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   dE = 0.001
      | 
      
   optical_matr_elem_model
      | Level: paragraph
      | 
      | Description: name of the *opticskp* model associated
      | 
      | 
      | Example: 
      |   optical_matr_elem_model = opticskp
      | 
      
   bias
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   bias = @V[0.0]
      | 
      
   load
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   load = @R[1e8]
      | 
      
   light_direction
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   light_direction = <vector_indicating_the_direction_of_light>
      | 
      
   light_intensity
      | Level: paragraph
      | 
      | Description: the vector which fixes the direction from where the light comes.
      | 
      | 
      | Example: 
      |   light_intensity = @x
      | 
      
   dye
      | Level: paragraph
      | 
      | Description: the kind of dye used in the cell
      | 
      | 
      | Example: 
      |   dye = N719
      | 
      
   mat
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   mat = TiO2mes
      | 
      
   porosity
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   porosity = 0.5
      | 
      
   TiO2
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   TiO2 = false
      | 
      
   generation
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   generation = dssc_generation
      | 
      
   values
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   values = ( 1000, 100, 1)
      | 
      
   poisson_simulation
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   poisson_simulation = dd
      | 
      
   doping_type
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   doping_type = donors
      | 
      
   doping_level
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   doping_level = 0.035
      | 
      
   force
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   force = (0.0, 0.0625,0.0)
      | 
      
   Young
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Young = 129
      | 
      
   Poisson
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Poisson = 0.349
      | 
      
   area_factor
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   area_factor = 0.1
      | 
      
   barrier_height
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   barrier_height = 3.0
      | 
      
   quadrature_rule
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   quadrature_rule = trapez
      | 
      
   flavour
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   flavour = relaxation
      | 
      
   monitor
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   monitor = true
      | 
      
   symmetry
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   symmetry = cylindrical
      | 
      
   ls_max_step
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   ls_max_step = 4
      | 
      
   discretization
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   discretization = fem
      | 
      
   integration_order
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   integration_order = 2
      | 
      
   Dirichlet_bc_everywhere
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   Dirichlet_bc_everywhere = true
      | 
      
   solution_method
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   solution_method = general
      | 
      
   solver
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   solver = krylovshur
      | 
      
   kp_model
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   kp_model = 8x8
      | 
      
   beta
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   beta = 1.0
      | 
      
   searchpath
      | Level: paragraph
      | 
      | Description: 
      | 
      | 
      | Example: 
      |   searchpath = ../../materials
      | 
   
   Modules/Device
      | Level: Super-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Device
      |    {
      |      meshfile = bulk.msh
      |    
      |      Region bulk 
      |      {
      |        material = Si
      |    
      |        Doping
      |        {
      |          Nd = 1e16
      |          type = donor
      |        }
      |      }
      |    }
      | 
      
   Modules/driftdiffusion
      | Level: Super-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Module driftdiffusion
      |    {               
      |      plot = (Ec, Ev, eQFermi, hQFermi, ContactCurrent)
      |    
      |      Contact anode { voltage = $Vb }
      |      Contact cathode { }
      |    }
      |
      
   Modules/sweep
      | Level: Super-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Module sweep
      |    {
      |      solve = driftdiffusion
      |      variable = $Vb
      |      start = 0.0
      |      stop = 1
      |      steps = 10
      |      plot_data = true
      |    }
      |
      
   Modules/Physics
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Physics
      |       {
      |         particle_density
      |         {
      |           statistics = boltzmann
      |         }
      |     
      |         mobility field_dependent
      |         {
      |           low_field_model = doping_dependent 
      |         }
      |         recombination (srh, auger) {}
      |       }
      |
      
   Modules/particle_density
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     particle_density
      |         {
      |           statistics = boltzmann
      |         }
      |
      
   Modules/mobility
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     mobility field_dependent
      |         {
      |           low_field_model = doping_dependent 
      |         }
      |
      
   Modules/Solver
      | Level: Super-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Solver
      |        {
      |         nonlinear_solver = tiber
      |         nonlin_max_it = 15
      |         nonlin_rel_tol = 1e-12
      |         nonlin_step_tol = 1e-5
      |        }
      |
      
   Modules/Cluster
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Cluster semiconductor
      |       {
      |        regions = -passivation
      |       }
      |
      
   Modules/NonlinearSolver
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     NonlinearSolver linesearch
      |        {
      |         relative_tolerance = 1e-15
      |         step_tolerance = 5e-3
      |         max_iterations = 15
      |
      
   Modules/LinearSolver
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     LinearSolver petsc
      |       {
      |        tolerance = 1e-7
      |        max_iterations = 10000
      |        #xmonitor = true
      |        pc = composite
      |       }
      |
      
   Modules/trap
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     trap fixed_charge
      |       { 
      |        region = GaN/SiN
      |        Nt = 2.74e13
      |       }
      |
      
   Modules/Boundary
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Boundary substrate
      |       {
      |        regions = bottom
      |        material = GaN
      |        y-growth-direction = (0,0,0,1)
      |        z-growth-direction = (1,0,-1,0)
      |        x-growth-direction = (-1,2,-1,0)
      |       }
      |
      
   Modules/BodyForcelattice_mismatch
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     BodyForcelattice_mismatch 
      |        {
      |         x-growth-direction = (-1,0,1,0) 
      |         y-growth-direction = (-1,2,-1,0) 
      |         z-growth-direction = (0,0,0,1)
      |         reference_material = AlGaN
      |         structure = wz
      |         x = 0.2 
      |        }
      |
      
   Modules/HeatSource
      | Level: Sub-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     HeatSource constant 
      |        {
      |         regions = HotSpot
      |         H = 1e10
      |        }
      | 
      
   Modules/efaschroedinger
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module efaschroedinger
      |       {
      |        particle = el
      |        poisson_model_name = dd
      |        strain_model_name = strain # macrostrain
      |        name = quantum_el
      |        regions = quantum
      |        plot = (EigenFunctions, EigenEnergy, EnergyLevels)
      |        Solver
      |          { number_of_eigenstates = 10 # 30 }
      |        Physics
      |          { model = conduction_band }
      | 
      
   Modules/quantumdispersion
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module quantumdispersion
      |       {
      |        simulation_name = dispersion1D_el
      |        regions = all
      |        quantum_simulation = quantum_el
      |        min_eigenvalue_number = 0
      |        max_eigenvalue_number = 5
      |        wedge = half
      |        k_space_dimension = 1
      |        k1 = (0, 0.1, 0)
      |        number_of_nodes = (10)
      |        output_format = grace
      |        plot = k-space_dispersion
      |       }
      |
      
   Modules/quantumdensity
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module quantumdensity
      |       {
      |        name = dens_el
      |        regions = quantum
      |        plot = quantum_density
      |        k-space_dimension = 0
      |        k-space_basis = true
      |        k1 = (0, 0, 0.1)
      |        quantum_simulation = quantum_el
      |        degeneracy = 2
      |        refine_fraction = 0.20
      |        relative_accuracy = 0.01
      |        refine_k_space = true
      |        output_density_in_k_space = true
      |        uniform_refinement = false
      |        mesh_order = FIRST
      |        first_state = 3
      |        initial_eigenstates_number = 10
      |        analytic = true
      |       }
      |
      
   Modules/opticskp
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module opticskp
      |      {
      |       name = optics
      |       regions = quantum
      |       plot = (optical_spectrum_k_0 )
      |       initial_state_model = quantum_el
      |       final_state_model = quantum_hl
      |       initial_eigenstates = (0, 9)
      |       final_eigenstates = (0, 15)
      |       polarization = (0, 0, 1)
      |       Emin = 2.8
      |       Emax = 3.6
      |       dE = 0.001
      |      }
      |
      
   Modules/opticalspectrum
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module opticalspectrum
      |       {
      |        k_space_dimension = 2
      |        k-space_basis = true
      |        k1 = (0, 0, 0.1)
      |        k2 = (0, 0.1, 0)
      |        refine_fraction = 0.30
      |        relative_accuracy = 0.01
      |        refine_k_space = true
      |        number_of_nodes = (2, 2)
      |        wedge = quarter
      |        plot = (optical_spectrum)
      |        optical_matr_elem_model = opticskp
      |        polarization = (0, 0, 1)
      |        Emin = 3.0
      |        Emax = 5.0
      |        dE = 0.001
      |       }
      |
      
   Modules/dssc
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     model dssc
      |       {
      |        options
      |          {
      |           simulation_name = <name_of_the_model>
      |           plot(Potential, Density, Current, ContactCurrents)
      |          }
      |       ...
      |       }
      
   Contact/anode
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact anode 
      |       {
      |         type = ohmic
      |         voltage = $Vb[0.0]
      |       }
      |
      
   Contact/cathode
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact cathode
      |       {
      |             type = ohmic
      |             voltage = 0.0
      |       }
      |
      
   Contact/gate
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact gate
      |      {
      |       regions = (gate1, gate2)
      |       type = schottky
      |       work_function = 1.5272
      |       voltage = @Vg
      |      }
      |
      
   Contact/Base
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact Base
      |      {
      |       type = clamp
      |      }
      |
      
   Contact/left
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact left
      |      {
      |       type = heat_reservoir
      |       temperature = 300
      |      }
      |
      
   Contact/right
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact right
      |      {
      |       type = heat_reservoir
      |       temperature = 300
      |      }
      |
      
   Contact/heat_reservoir
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact heat_reservoir
      |      {
      |       region = substrate  
      |       temperature = 300
      |      }
      |
      
   Contact/substrate
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact substrate
      |       {
      |        type = surface_resistance
      |        r_surf = 0. 05 
      |        temperature = 300
      |       }
      |
      
   Contact/reservoir
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact reservoir
      |       {
      |        temperature = 300
      |       }
      |
      
   Contact/surface_resistance
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact surface_resistance
      |       {
      |        r_surf = 0.01
      |        temperature = 300
      |       }
      |
      
   Contact/Top
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact Top
      |       {
      |        type = surface_force
      |        force = (0.0, 0.0625,0.0)
      |       }
      |
      
   Contact/backcontact
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact backcontact
      |       {
      |        voltage = 0.0
      |        area_factor = 0.1
      |       }
      |
      
   Contact/dissipator
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact dissipator
      |       {
      |        r_surf = 0.5
      |        type = surface_resistance
      |        temperature = 300
      |       }
      |
      
      
   Modules/dssc_generation
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     model dssc_generation
      |       {
      |        options
      |          {
      |           regions = (<TiO2_region_name_1>, <TiO2_region_name_2>, ...)
      |           plot = (Distance, Generation)
      |           light_direction = <vector_indicating_the_direction_of_light>
      |          (example, illumination from x direction: light_direction = (1, 0, 0))
      |           light_intensity = @x
      |           dye = N719
      |          }
      |       ...
      |       }
      |
      
   Modules/sweep_gen
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     sweep_gen
      |       {
      |        simulation = (dssc_generation, dssc)
      |        variable = x
      |        values = (0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1)
      |        plot_data = true
      |       }
      |
      
   Modules/BodyForceconverse_piezo
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     BodyForceconverse_piezo
      |       {
      |        poisson_simulation = dd
      |       }
      |
      
   Modules/elasticity
      | Level: Super-Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module elasticity 
      |       {
      |        Physics 
      |          {
      |           Stiffness isotropic
      |             {
      |              Young = 129 Poisson = 0.349
      |             }
      |           Contact Base
      |             {
      |              type = clamp
      |             }
      |           Contact Top
      |             {
      |              type = surface_force force = (0.0, 0.0625,0.0)
      |             }
      |          }
      |       }
      |
      
   Modules/Stiffness
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Stiffness isotropic
      |       {
      |        Young = 129 
      |        Poisson = 0.349
      |       }
      |
      
   Modules/Selfconsistent
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Module Selfconsistent
      |       {
      |        flavour = relaxation
      |        simulations = (dd,tt)
      |        abs_tolerance = 1e-5
      |        monitor = true
      |       }
      |
     
   Region/bulk
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region bulk 
      |       {
      |         material = Si
      |     
      |         Doping
      |         {
      |           Nd = 1e16
      |           type = donor
      |         }
      |       }
      |
      
   Region/n_side
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region n_side 
      |      {
      |        Doping
      |        {
      |         Nd = 1e18
      |         type = donor
      |        }
      |      }
      |
      
   Region/p_side
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region p_side 
      |       {
      |         Doping
      |         {
      |          Nd = 1e18
      |          type = acceptor
      |         }
      |       }
      |
      
   Region/barrier
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region barrier
      |       {
      |         material = AlGaN
      |         x = 0.14
      |        }
      |
      
   Region/barrier_doped
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region barrier_doped
      |       {
      |        mesh_regions = (barrier_dop_s, barrier_dop_d)
      |        material = AlGaN
      |        x = 0.14
      |        Doping 
      |          {
      |           density = 1e21
      |           type = donor
      |           level = 0.026
      |          }
      |       }
      
   Region/buffer_doped
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region buffer_doped
      |       {
      |        mesh_regions = (buffer_dop_s, buffer_dop_d)
      |        Doping 
      |          {
      |           density = 1e21
      |           type = donor
      |           level = 0.026
      |          }
      |       }
      |
      
   Region/cap
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region cap
      |       {
      |        Doping 
      |          {
      |           density = 5e18
      |           type = donor
      |           level = 0.026
      |          }
      |       }
      |
      
   Region/cap_doped
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region cap_doped
      |       {
      |        mesh_regions = (cap_dop_s, cap_dop_d)
      |        Doping 
      |          {
      |           density = 1e21
      |           type = donor
      |           level = 0.026
      |          }
      |       }
      |
      
   Region/passivation
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region passivation
      |       {
      |        mesh_regions = (passiv1, passiv2, passiv3, passiv4)
      |        material = SiN
      |       }
      |
      
   Region/Well
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region Well
      |       {
      |        material = GaN
      |       }
      |
      
   Region/porous
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region <name of the porous region>
      |       {
      |        mat = TiO2mes
      |        porosity = 0.5
      |       }
      |
      
   Region/electrolyte
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region <name of the electrolyte region>
      |       {
      |        mat = TiO2mes
      |        TiO2 = false
      |       }
      |
      
   Region/Nanowire 
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region Nanowire 
      |       {
      |        doping_type = donors
      |        doping_level = 0.035
      |       }
      |
      
   Region/Air 
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region Air
      |       {
      |        material = Air
      |       }
      |
      
   Region/oxide
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region oxide
      |       {
      |        material = SiO2
      |       }
      |
      
   Region/collector
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Region coll # the collector region
      |       {
      |        type = donor
      |        Nd = 5e15
      |       }
      |
   
   Contact/anode
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact anode 
      |       {
      |         type = ohmic
      |         voltage = $Vb[0.0]
      |       }
      |
      
   Contact/cathode
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact cathode
      |       {
      |             type = ohmic
      |             voltage = 0.0
      |       }
      |
      
   Contact/gate
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact gate
      |      {
      |       regions = (gate1, gate2)
      |       type = schottky
      |       work_function = 1.5272
      |       voltage = @Vg
      |      }
      |
      
   Contact/Base
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact Base
      |      {
      |       type = clamp
      |      }
      |
      
   Contact/left
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact left
      |      {
      |       type = heat_reservoir
      |       temperature = 300
      |      }
      |
      
   Contact/right
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact right
      |      {
      |       type = heat_reservoir
      |       temperature = 300
      |      }
      |
      
   Contact/heat_reservoir
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |    Contact heat_reservoir
      |      {
      |       region = substrate  
      |       temperature = 300
      |      }
      |
      
   Contact/substrate
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact substrate
      |       {
      |        type = surface_resistance
      |        r_surf = 0. 05 
      |        temperature = 300
      |       }
      |
      
   Contact/reservoir
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact reservoir
      |       {
      |        temperature = 300
      |       }
      |
      
   Contact/surface_resistance
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact surface_resistance
      |       {
      |        r_surf = 0.01
      |        temperature = 300
      |       }
      |
      
   Contact/Top
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact Top
      |       {
      |        type = surface_force
      |        force = (0.0, 0.0625,0.0)
      |       }
      |
      
   Contact/backcontact
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact backcontact
      |       {
      |        voltage = 0.0
      |        area_factor = 0.1
      |       }
      |
      
   Contact/dissipator
      | Level: Block
      | 
      | Description: 
      | 
      | 
      | Example: 
      |     Contact dissipator
      |       {
      |        r_surf = 0.5
      |        type = surface_resistance
      |        temperature = 300
      |       }
      |
      
