# tiberCAD TODO list

This is the TODO list for tiberCAD software project. Format is following
https://github.com/todomd/todo.md
Here we put everything to be done in the tiberCAD core, while each module can have its
own TODO list. For historical reason, some module-specific TODO's are left in this list, however.

### Todo

- [ ] reorganize Makefiles and add out-of-src build
- [ ] reorganize and unpack stuff in Utils.h
- [ ] eliminate tiberCAD.in in bin, transfer to TiberCAD.sh.in and link  
- [ ] eliminate interlayer for libMesh system and solver classes (and solver re-implementations)
- [ ] add more useful common API for projecting quantum states on bases/orbitals
- [ ] change input parser so that it reads first into a string (with all includes) and then parses 
- [ ] change get_solution to have source argument, to allow different meshes more easily  
- [ ] allow reading individual meshes in modules
- [ ] get rid of read_database(), put database access into get_option() and get_parameter()  
- [ ] throw exception if get_option() gets a string with "$"  
- [ ] in 1D/2D check if mesh is along x or in x/y  
- [ ] improve/document tunneling current (analytical model)  
- [ ] implement inspect output  
- [ ] in DriftDiffusion, put calculated currents into contacts  
- [ ] implement model for Zb valence band mass at small strain  
- [ ] allow for more complex aomistic structures, like wz-zb

### In Progress

- [ ] eliminate Tensor1 and Tensor2 in favor of libmesh classes
- [ ] add possibility to have independent meshes  
- [ ] add transient solver 


### Done ✓


- [x] eliminiate RotatedCrystal (transfer to BulkCrystal)  
- [x] get rid of all the bad static pointers to 'this' in modules  
- [x] control/change get_solution() for CELL data (return data on each point)  
- [x] add something like @define, @if etc. in InputParser  
- [x] Complete the Elasticity module  
- [x] add possibility to read RealGradient type from input file  
- [x] add the \*operator between tensor.h and tensor_value.h related types.  
- [x] put plot keyword into modules  
- [x] put mesh, dimension into $Device  
- [x] review approach for variables  
- [x] implement surface recombination  
- [x] implement more trap types  
- [x] cutoff in embracing region  
- [x] implement doping profile/external doping  
- [x] implement Fowler-Nordheim for contacts  
- [x] improve Sweep to make 'logarithmic' scale sweeps, change meaning of min/max steps  
- [x] change electrical contact for the case of two or more materials touching the contact  
- [x] change surface charge model so it can work at the interior  
- [x] improve drift-diffusion for electron-only calculations  
- [x] immplement quantum Hamiltonian model for X and L valleys  
- [x] put field dependent mobility into jacobian  
- [x] throw exception when vb DOS mass negative  
