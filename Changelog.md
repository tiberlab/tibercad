# Changelog

All notable changes to tiberCAD should be reported here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Modules can have their own changelog file.


## [Unreleased]

## [v3.5.0] - 2024-12-22

### Added

- implemented reordering of states when calculating dispersions based on projection on nearby k-points
- added new module pvmodule  
- added new module wateringress  

### Fixed

- moved efaschroedinger and opticskp to proper modules  

### Changed

- changed several filenames in preparation of migrating to github  


# [v3.4.0] - 2024-10-03

### Added
- upgraded to libmesh 1.7.1, PETSc/SLEPc 3.17, cuda 11.5  
- several fixes  
- can now provide crystal directions as planes (hkl) or directions
	  [utw], both Miller or Miller-Bravais. Also, only one or two
	  directions need to be given  
- can now specify values when to plot data in Sweep  


## [3.3.0] - 2024-04-17

- implemented Tan scheme in uptight  
- ETB parameters are organized in subdirectories for different schemes  
- TMM has now internal source, can be coupled to DD  
- NEGF module can use the new quasi-equilibrium density calculation  

## [3.2.0]


### Added
- upgraded to libmesh 1.6.3, PETSc/SLEPc 3.16.1, boost 1.77, cuda 10.2  
- several fixes  
- added some features for k-space and crystal symmetries  
- new module for TMM  
- added interfacing to SLEPc in ETB module (solver = slepc)  
- intraband terms now included in kp optical calculations  
- in kp, can now project solutions onto basis states to get
	  composition  

## [3.1.0]


### Added
- can now read GMSH msh format version 4.1  
- new data import module, can read 1/2/3D ASCII data  
- can now handle 2D materials for atomistic simulations  
- improvements for k space generation (dispersion and integration)  
- implemented unfolding of supercell dispersions (for wurtzite and
	  simple cubic/tetragonal currently)  
- empirical tight binding works now for any orbital combinations  
- new examples for tight binding, and for doping profiles / alloy
	  profiles  
- new module to create quantum espresso input files  
- alloy statistics extraction can now include APT detection efficiencies
	  and use binning based on counting atoms  
- mp-DD functional and mostly complete  
- nonlinear piezoelectricity implemented for zb and wz (PRB 84 and 88)  


## [3.0.2]


### Fixed
- fixed a bug in dssc_generation for cylindrical symmetry  


## [3.0.0]


### Added
- upgraded to libmesh 1.0.0, PETSc 3.6.2, SLEPc 3.6.2  
- improved MPI parallelization: can assign communicators
	  to devices and modules  
- k-integration can parallelize automatically  
- contacts can be defined embedded in single regions  

### Fixed
- fixed a bug in dssc_generation  
- fixed a bug for strain corrected bulk band parameters  


## [2.6.0]


### Added
- kp based band parameters now explicitly uses all available bands to
	  calculate the carrier densities, instead of extracting effective
	  parameters for a single parabolic band (which has been done only for
	  holes, anyway)  
- quaternaries can now be used for atomistic simulations  
- optics modules can now use Gaussian line shape (line_shape =
	  gaussian | lorentzian)  
	

## [2.5.0]


### Added
- complete reorganization of handling of the band parameters in drift
	  diffusion. Allows now for different DOS implemented as modules  
- implemented several models for the simulation of organic materials,
	  e.g. Gaussian DOS, Langevin recombination  
- Valence force field is included  
- Empirical tight binding is included  

### Fixed

- missing factor of 2/3 added in spontaneous emission  

## [2.3.0]


### Added

- optics now calculates also stimulated emission and gain spectra  
- recombination from optics can be included selfconsistently in
	  drift-diffusion  
- EFA can now treat electrons and holes in a unified way  
- k-space is measured in 1/nm  
- can now define a default for doping in Device block  
- can use Projection block in modules to project solutions on a
	  regular tensor grid  
- atomistic structure can calculate alloy statistics  
- trap occupation is now based on both quasi Fermi levels  
- can have generation rates between SRH trap states and bands  
- output for recombination rates is given for electrons and holes
	  distinctly  
- can have variable alloy fractions, based on external input or with
	  uniform random distribution (plot_alloy_composition = true for
	  plotting)  

### Fixed

- interface models for drift-diffusion now correctly consider the
	  min/max of the band edges and densities from the two sides  
- extended boundary condition in elasticity module works correctly  
- bug fix in TAT model  


## [2.2.0]


### Added

- logfile is now written in 'resultpath' directory  
- bug fix for recombination models  
- new option in sweep 'ignore_failure = true' will ignore
	  solver failures and simply exit the sweep  
- integer and double vectors accept ranges as a-b, or a:b, or a:s:b,
	  e,g, (0, 0.1, 0.2:0.1:0.6, 1-3)  
- print some info on errors and warnings at the end  
- new flag -i for interactive mode, which will stop at errors
	  and warnings  
- input file macros @define, @ifdef, @ifndef, @else, @endif  
- in driftdiffusion can now apply optical generation at interface  
- bug fixes in elasticity module, shape deformation is now working correctly  
- implemented local version of BBT and TAT (Hurkx model)  
- implemented local version of Schottky barrier tunneling  
- a bug fix in optics module, and can now specify the effective
	  refractive index in the input file  
- Quaternary alloys of type (A_yBD)_{x}(C_zBD)_{1-x} are now supported  
- write_boundary_mesh=true in device section will write the region
	  boundaries as VTK file (for dim > 1)  


## [2.1.0]


### Added

- improved automatic guess for eigenenergies in quantum
	  simulations  
- new submodels conduction_band and valence_band in driftdiffusion,
	  allowing for better control of band parameters (or band_parameters
	  to control them usaing a single block)  
- extended syntax for plot variables: can now switch off variables by
	  writing plot = (foo, -bar), e.g. plot = (all, -Strain)  
- can now use single band approximation for holes  
- can now define interfaces by specifying material or region pairs (or
	  a combination of both) : mat1/mat2, region1/mat2. External boundary
	  is specified by a dash '-': Si/-  
- improved generation model for DSC based on Lambert-Beer approximation  
- improved boundary conditions for DSC including Butler-Volmer for
	  ionic species and kinetic rate for electrons  
- Multi trapping model to include in an effetive way the exponential band
	  of traps below TiO2 conduction band edge  

	

## [2.0.0]


### Added

- new input file syntax  
- new physical models (polarization, particle_density)  
- new modules (elasticity)  
- simplified input file syntax for EFA module (incorporates now
	  quantum density)  
- improved selfconsistency Quantum/Drift-Diffusion based on
	  corrector/predictor  
- windows version can open file dialog to open input file  
- can now selectively exclude regions by saying e.g. regions = (all, -bulk)  

### Changed

- xx_rel_tol -> relative_tolerance  
- xx_abs_tol -> absolute_tolerance  
- xx_max_it  -> max_iterations  
- pc_type -> preconditioner  
- ksp_type -> method  

## [1.X.X]

### Added

- added new materials  
- simple optical generation model  
- New 'Options' block in $Device section of input file for
	  specification of common options for device regions  
- save and load of solutions in drift-diffusion  



### Fixed

- fixed error in database for GaN, InN direct recombination  
- fixed bug regarding variables in alloys  
- fixed sign issue for zincblende deformation potentials  
- fixed bug in mesh reading (now  it's fast)



## [1.2.2]


### Added

- Auger recombination  
- Surface recombination for interface boundary model  
- self consistency Schroedinger-Poisson/Drift Diffusion, based on 
predictor-corrector scheme with embracing region  
- physical models can be restricted to a subset of the simulation domain  
- string labels can be used in GMSH as physical region identifiers instead of numbers  
- more meaningful error messages from input parser  

### Fixed

- fixed missing error messages in Windows version  
- fixed bug in field dependent mobility model  
- fixed deformation potentials for nitride semiconductors  


