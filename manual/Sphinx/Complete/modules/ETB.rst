..   <marker>

.. _ETBTheory:


Empirical Tight Binding
=================================



In  tiberCAD,  it  is  possible  to  perform atomistic quantum  calculations in  the  framework  of  Empirical Tight Binding (ETB):  eigenstates, eigenfunctions and  quantum  density of  a  given system   can  be  obtained  by solving a tight-binding Hamiltonian by means  of the module:

* **Module**  ``empirical_tb``

The  optical properties  are  calculated by  the  module 

* **Module**  ``opticstb``



A typical example of tight-bind calculations is the following  ::

  Module empirical_tb
  {

    regions = Quantum
    name = tb
    atomistic_structure = tb
    Harrison_scaling = true 

    strain_simulation = strain
    potential_simulation =  driftdiffusion

    sparse_format = full
 
    plot = (tbstates, ProbabilityDensity)
    jmol_output_format = cube

  

     Solver
      {
        #load_states = true 
        #load_path = output3
         num_valence_eigenvalues = 2 
        num_conduction_eigenvalues = 2 
        long_tolerance = 1e-6
        guess_valence = -1.4 
      }

  }#empirical_tb




Atomistic  generator
-----------------------

In the prevoius  example, TB  calculations  are applied  to  the  physical  region ``Quantum`` (usuallly a  Cluster) ::

  regions = Quantum   


An  atomistic structure must be  coupled with  this physical region, in  order  to  apply the TB  model ::

  atomistic_structure = tb

This  is  made through the **Atomistic  generator**   a  quite  powerful tool  of  tibercad which  allows  to  generate an atomic  basis associated  with the finite element mesh which belongs to  a  given  physical region, based  on  the material specifications and the growth directions defined  for  that region.
 
For  example, the  atomistic  structure ``tb``  is created  by  the  **Atomistic  generator** through the  following block of instructions contained in  the  *Device*  section  ::

  Atomistic tb
	{
	  reference_region = barrier_left
	  regions = Quantum 
	  passivation = yes
	  print = (xyb, xyz, gen, tgn)
	}


The internal atomistic generator has the following features

* It supports any crystal structure with fcc, bcc, cubic and hexagonal Bravais lattice.
  Informations about basis and primitive vectors is taken from the database, and new materials can be added without changing the code.


* It performs hydrogen passivation for any crystal structure supported, 



The  definition of  a  *reference region* is  needed. It provides the  specifications on the material and parameters (growth direction, molar fraction, crystal phase) we will use to build the crystal. 

Note that we assume a pseudomorphic growth, i.e. when  (which is  usyually the  case) we generate  a heterostructure, this has the same lattice constants of the  reference region.

Taht is, only the atomic species are changed according to the different materials specifications, while the  lattice is  the reference material lattice. 
 
The structure may  be eventually relaxed by  projecting strain  calculated with a continuum elasticity model, or by  minimizing  the  energy  of  the  system applyin  a  valence force field (VFF) model.


Options
^^^^^^^^^^^^^^^^^^^^^^^

The  following kewords may  be  defined in *Atomistic*  block


 ``regions`` : string
    defines the phisical regions or  Cluster where to  build  the  atomistic  structure

 ``reference_region`` : string 
    defines the name of the region whose  material  and  growth directions  are  taken  as  a  reference for  lattice  construction

 ``reference_material`` : string 
    defines a general reference  material  and  growth directions that are taken as a  reference for  lattice  construction

 ``passivation`` : boolean
    if true an hydrogen passivation is  applied to the structure. 
    Hydrogen atoms are inserted at a defined distance in the same direction of the original crystal bond.

 ``hydrogen_distance`` : real
    defines the bond length at which Hydrogens are placed from the surface atoms (default 1.2 Angstrom)

 ``periodic`` : boolean
    if **true** makes the structure periodic. For instance for a 3D mesh the structure is made periodic in all 3 cartesian directions.
    Limitations: the geometry must be a simple box with sides oriented along the cartesians directions.

 ``translation`` : real vector
    defines a translation vector as (dx, dy, dz) used to move the atomistic structure w.r.t. the mesh regions. 
    It is useful for fine alignments of atoms to the mesh. 
 
 ``preserve`` : string
    **Not active yet**. In principle should be used to ensure that all atoms of unit cell or conventional cell are kept in the structure

 ``load or load_structure`` : string
    specifies a (relative) path/file for loading of an external structure. valid formats are tgn, xyz, gen.
    The structure is trimmed out the regions on which the atomistic structure is defined. 
    Note: when loading xyz files, the supercell vectors must be defined in the input file.

  ``lattice_vectors`` : list of floating
    specifies the periodicity  vectors used  in  a  supercell  calculation.
    The format is (x1,y1,z1,x2,y2,z2, x3,y3,z3),  values  are  in  A, e.g.( 52.728 , 0.0, 0.0 , 0.0, 45.664 , 0.0, 0,0, 42.758 )

 ``y_length`` : real
    defines the lateral supercell size (in Angstrom) in 1D (quantu wells) and 2D (nanowires) structures.

 ``z_length`` : real
    defines the lateral supercell size (in Angstrom) in 1D and 2D structures.

 ``random_alloy`` : boolean
    defines generation of a random alloy structure (default = **false**). 

 ``random_generator_seed`` : integer
    an integer defining the seed. Useful to produce exactly the same random structure.   

 ``fix_mean_alloy_concentration`` : boolean
    if **true** the number of substituted atoms is as close as possible to the nominal molar fraction of the alloy.
    if **false** each atom is substituted with a probability proportional to the molar fraction. 

 ``clustering`` : boolean
    if **true** an alloy with clustering is generated. This is done by placing a fraction of seed atoms (cluster_seed) with uniform 
    distribution, then the remaing fraction of atoms is placed biased with a larger probability near existing ions.  

 ``cluster_seeds`` : real
    defines the fraction of ions placed with uniform distribution. Must be >0. (default = 0.02). 

 ``extract_alloy_statistics`` : boolean
    if **true** the local concentration is computed by averaging over small local spheres of radius defined by control_volume_radius.

 ``control_volume_radius`` : real
    defines the radius of a sphere on which local concentration is computed (default = 0.5 nm).     
   
 ``reference_atom`` : integer
    sets the atom label (typically 1=cations, 2=anions) where local concentration is computed     

 ``plot_alloy_statistics`` : boolean 
    if **true**, the local concentration is mapped on the mesh and plotted in a mesh format.

 ``meshdata_format`` : string
    defines the output format of the alloy statistics (default = "vtk")

 ``print`` : string
   *xyz*:  xyz file  format;
   *xyb*: modified xyz format with additional columns containing the bond map;
   *gen*: gen file  format;
   *tgn*: tibercad format which combines informations on the mesh and the atoms.

Module empirical_tb 
-----------------------


Module options
^^^^^^^^^^^^^^^^^^^^^^^


The following options influence the behaviour of the **Module**  ``empirical_tb``:

 ``atomistic_structure`` : string  
    defines the name of the atomistic structure which  will  be  used  by ETB  to  create the  TB Hamuiltonian. It  must be  the name  of  one  structure created  through an **Atomistic** block in  **Device**  section. 
   

 ``potential_simulation`` : string
    defines the name of the simulation (e.g.  ``driftdiffusion``) that can provide electric potential

 ``strain_simulation`` : string 
    defines the name of the simulation (e.g. ``elasticity``) that can provide elastic strain

 ``regions`` : string 
    defines the regions associated to this ETB simulation

 ``Harrison_scaling`` : boolean
    if  **true** , scaling  of  ETB  parameters is  applied (usually in  presence  of  material  deformation which  causes atom  displacement  from  equilibrium  position).
    This  option is  by deafult  set to **true**   whenever  a  strain  simulation is  performed on  the  system, otherwise is  by  default  **false**.



 ``sparse_format`` : string 
    Possible values are ``full`` and  ??



Solver section
^^^^^^^^^^^^^^^^^^^^^^^


The Solver section of the **Module**  ``empirical_tb`` contains the following options:

 ``num_valence_eigenvalues`` : integer 
    defines the number of eigenvalues and eigenfunctions to be found.

 ``num_conduction_eigenvalues`` : integer 
    defines the number of eigenvalues and eigenfunctions to be found.

 ``load_states`` : boolean 
    if **true**, calculated states are  loaded from  the  file specified by  ``load_path``

 ``load_path`` : string 
    path of  the file  from  which  tb  states  can  be loaded

 ``long_tolerance`` : double 
    tolerance...




Output
^^^^^^^^^^^^^^^^^^^^^^^

The available output variables for **Module**  ``empirical_tb`` , to be specified in the plot option, are the following:

 ``tbstates`` :
   states eigenfuctions in  *cube*  format

 ``MeshStatesNodes``: 
    plots  staes  on  the  mesh. MeshStatesNodes is not supported in 1D calculations 

 ``MeshStates``: 
    plots  states  on  the  mesh in 1D calculations (ee.g. Quantum wells).  



By  deafault  the list of  eigenvalues[eV] is printed  in the output  directory 


Module opticstb
----------------------

The **Module** ``opticstb`` implements the  calculation of optical properties based on  the  tb states; in general, the optical  matrix elements are calculated from the tb models specified in ::

  initial_state_model = tb2
  final_state_model = tb2


The main parameters are:

 ``initial_state_model``:
     tb simulation associated to the initial state of the optical transition

 ``final_state_model`` :
    tb simulation associated to the final  state of optical transition

 ``Emin``:
    minimum  energy 

 ``Emax``:
    max  energy

 ``dE``:
    energy step



 ``plot`` :
    *optical_spectrum_k_0* to  select  spectrum in *k=0* (:math:`\Gamma` point) 
    *matrix_elements* to  select calculation of  the ............


Integrated spectrum
^^^^^^^^^^^^^^^^^^^^^^^


.....................




Example
-----------------------


In this example we will see how to calculate quantum states of a GaN **quantum dot**  with a tight-binding model based on a sp3s*d parameterizaation. Strain  due  to  the  lattice  mismatch is taken  in  account  with  elasticity  calculations  and  the  related  atom  displacements are projected  on  the  atomic  structure.
Potential  profiles including  piezo and  pyro  polarizations  are  calculated  by  solving  Poisson equation and  projected  as  well  on  the  TB Hamiltonian.

First,  the  Device  structure  is  defined ::

  Device qdot
  {
   meshfile = box_in_cube.msh 
   mesh_units = 1e-9
   dimension = 3

    y-growth-direction = (-1,0,1,0)
    x-growth-direction = (-1,2,-1,0)   
    z-growth-direction = (0,0,0,1)

  Region ball
   {
    material = GaN  

    Doping {
      Nd = 1e15
      type = donor
      Ed = 0.025
    }

   }

  Region  qbox
   {
    material = AlN  
 
    Doping {
      Nd = 1e15
      type = donor
      Ed = 0.025
    }
    
  }


  Region  intrinsic
  {
    material =  AlN  
   
    Doping {
      Nd = 1e16
      type =  donor
      Ed = 0.025
    }
    
  }

  Region  nside
  {
    material =  AlN 
   
    Doping {
      Nd = 1e18
      type =  donor
      Ed = 0.025
    }
    
  }


  Region  pside
  {
    material =  AlN  
   
    Doping {
      Nd = 1e18
      type =  donor
      Ed = 0.025
    }
    
  }


  Cluster atomistic
  {
    regions=(ball, qbox)
  }


  Atomistic tb1
  {
    regions = atomistic 
    reference_region = nside  # barrier1 
    #path = tb1.gen
    passivation = yes
    #preserve = conventional
  }

}



A   GaN  quantum  dot (QD) is  surrounded by  an AlN  quantum  box:   these two  regions,  named ball  and  qbox, constitute  the atomistic cluster.  Based  on  the  finite  element grid defined  by  this cluster,  an  atomistic  structure  is  build,  as  defined  in **Atomistic**  block  of  the  Device  section.
The  reference  lattice  used  to  build  the  structure is  defined by ::

  reference_region = nside 
  

to  be that of the   AlN  material of  one  of  the other  regions  which  form the  Device: two  doped and  one  intrinsic  AlN  regions.

Strain  due  to  lattice  mismatch  between GaN  and AlN (taken  as  reference  substrate) is  taken  into  account  through the Module  elasicity ::

  Module elasticity
  {

  name = str  
  regions = all

  plot = (Strain, Stress, Displacement )

  non_linear_strain = true
  shape_error = 0.1


  Physics
  {
    body_force lattice_mismatch
    {
      reference_material = AlN
     
      y-growth-direction = (-1,0,1,0)
      x-growth-direction = (-1,2,-1,0)   
      z-growth-direction = (0,0,0,1)

    }

  }

  Contact cathode
  {type = clamp}

  }



Tight  binding  calculation  of  quantum  states  in  the  system  made  by  the  QD  and  the  AlN  surrounding  box is  defined  by  the  Module empirical_tb ::

  Module empirical_tb
  {
   regions = atomistic 
    name = tb

    atomistic_structure = tb1
    potential_simulation = dd
    strain_simulation = str  
    sparse_format = full	
    plot = (tbstates)

    Solver
     {
      num_valence_eigenvalues = 1 
      num_conduction_eigenvalues = 1
      long_tolerance = 1e-4  
      }

  }


Note  that  the  results of potential (from driftdiffusion)  and  strain (from  elasticity) simulations   are  recalled,  to  project  the  correct  potential  profile  and  atom  displacement  to  the  TB  Hamiltonian.
In  output,  the ground  states  of  conduction  and  valence  bands  are saved  together  with  the  wavefunctions.

Based  on  these  results, *optical  properties*  are  calculated through  the Module  opticstb ::

  Module opticstb
  {
   name = opt
   regions =  atomistic 

   initial_state_model = tb2
   final_state_model = tb2
  
   compute_strengths = true

   plot = (matrix_elements, optical_spectrum_k_0)
   output_format = grace

   Emin = 4.25
   Emax = 4.75
   dE = 0.001

  }



The  initial  and  final  state  for  the  spontaneous  emission  recombination are  defined  by ::

  initial_state_model = tb2
  final_state_model = tb2

where ''tb2'' is  the  name  of  the  empirical tb  simulation  which  performs  calculations 
Usually,  this  is  the  same  for  electrons and holes.

The  energy range  for  the  calculation of  the  optical  spectrum  is  defined  through ::

  Emin = 4.25
  Emax = 4.75
  dE = 0.001



 
