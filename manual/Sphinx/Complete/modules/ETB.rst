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

This  is  made through the **Atomistic  generator**,   a   tool  of  tibercad which  allows  to  generate an atomic  basis associated  with the finite element mesh which belongs to  a  given  physical region, based  on  the material specifications and the growth directions defined  for  that region.
 
For  example, the  atomistic  structure ``tb``  is created  by  the  **Atomistic  generator** through the  following block of instructions contained in  the  *Device*  section  ::

  Atomistic tb
	{
	  reference_region = barrier_left
	  regions = Quantum 
	  passivation = yes
	  print = (xyz, gen, tgn)
	}


The internal atomistic generator has the following features

* It supports any crystal structure with fcc, bcc, cubic and hexagonal Bravais lattice.
  Information about basis and primitive vectors is taken from the material database.

* It performs hydrogen passivation for any crystal structure supported. 

* For the  generation  of heterostructures  a pseudomorphic growth is assumed, with the same lattice constants of the  defined reference region.



The  definition of  a  *reference region* is  needed. It provides the  specifications on the material and parameters (growth direction, molar fraction, crystal phase) which  will be used to build the crystal structure. 

In this way,  only the atomic species are changed according to the different materials specifications, while the  lattice is  the reference material lattice. 
 
The structure may  be eventually relaxed by  projecting strain  calculated with a continuum elasticity model, or by  minimizing  the  energy  of  the  system applying  a  valence force field (VFF) model.


Options
^^^^^^^^^^^^^^^^^^^^^^^

The  following kewords may  be  defined in *Atomistic*  block


 ``regions`` : string
    defines the Phisical Regions or  Cluster where to  build  the  atomistic  structure

 ``reference_region`` : string 
    defines the name of the region whose  material  and  growth directions  are  taken  as  a  reference for  lattice  construction

 ``reference_material`` : string 
    defines a general reference  material  and  growth directions that are taken as a  reference for  lattice  construction (may be useful if one chooses a  reference material not present in  any device region)

 ``passivation`` : boolean
    if true an hydrogen passivation is  applied to the structure. 
    Hydrogen atoms are inserted at a defined distance in the same direction of the original crystal bond.

 ``hydrogen_distance`` : real
    when *passivation = true* , defines the bond length at which Hydrogen atoms  are placed from the surface atoms (default 1.2 Angstrom

 ``periodicity`` : list of boolean
    add periodicity along one or more  of the three  directions *x,y,z*, depending on  the mesh dimension: for a 1D mesh, add periodicity in *x* direction (*periodicity = true*); for a 2D mesh, add periodicity along *x* and/or *y*  directions (e.g. *periodicity = (true,true)*; for a 3D mesh, add periodicity along *x,y* and *z* (e.g. *periodicity = (true,true,true)*. More examples in  the  following.

 ``translation`` : real vector
    defines a translation vector as *(dx, dy, dz)* used to move the atomistic structure w.r.t. the mesh regions. 
    It is useful for fine alignments of atoms to the mesh. 
 
 ``load_structure`` : string
    specifies a (relative) path/file for the loading of an external structure. Valid formats are *tgn, xyz, gen*.
    The external structure is trimmed outside  the regions on which the atomistic structure is defined in the input file. 
    Note: when loading *xyz* files, the supercell vectors must be defined  through the keyword *lattice_vectors* (see below)

 ``lattice_vectors`` : list of floating
    specifies the periodicity  vectors used when importing an external atomic structure in *xyz* format (see *load_structure*). The format is (x1,y1,z1,x2,y2,z2, x3,y3,z3),  values  are  in  A, e.g.( 52.728 , 0.0, 0.0 , 0.0, 45.664 , 0.0, 0,0, 42.758 ). 

 ``supercell_size_y`` : real
    defines the lateral supercell size on *y-axis* (in Angstrom) in 1D (quantum wells) and 2D (nanowires) structures.

 ``supercell_size_z`` : real
    defines the lateral supercell size on *z-axis* (in Angstrom) in 1D and 2D structures.

 ``random_alloy`` : boolean
    if **true**,  a random alloy structure is generated (default = **false**).More examples in the following. 

 ``random_generator_seed`` : integer
    an integer defining the random seed. Useful to produce exactly the same random alloy structure in a  set  of  executions.   

 ``fix_mean_alloy_concentration`` : boolean
    if **true** the number of substituted atoms in a random alloy is as close as possible to the nominal molar fraction of the alloy.
    If **false** each atom is substituted with a probability proportional to the defined nominal molar fraction. The resulting actual number of substituted atoms may differ from the nominal  concentration value.

 ``clustering`` : specie name 
    Alloy models with clustering on the given atomic specie are produces. This is done by placing a fraction of seed 
    atoms (cluster_seed) with uniform distribution, then the remaing fraction of atoms is placed biased with a larger 
    probability near existing ions.  

 ``cluster_seeds`` : real
    defines the fraction of ions placed with uniform distribution. Must be > 0. (*default = 0.02*). 

 ``control_volume_radius`` : real
    defines the radius of a sphere on which local concentration is computed (default = 0.5 nm).     
   
 ``reference_atom`` : integer
    sets the atom label (typically 1=cations, 2=anions) where local concentration is computed. Default is cations.    

 ``meshdata_format`` : string
    defines the output format of the alloy statistics (default = "vtk") plotted when *plot_alloy_statistics = true*.

 ``print`` : string
    print the  generated atomic  structure with  one  of  the  following formats: *xyz*:  xyz file  format; *gen*: gen file  format; *tgn*: tibercad format which combines informations on the mesh and the atoms.

 ``alloy_statistics`` : block 
    if defined, the local concentration is computed by averaging over small local spheres of radius 
    defined by control_volume_radius. Example

  alloy_statistics
  {
    control_volume_radius = 1.0 
    plot_alloy_composition = true
  }
 
 ``control_volume_radius`` : real
   sets the radius of the control volume. 

 ``plot_alloy_composition`` : boolean
   if **true** the alloy composition is mapped on the mesh and plotted in a vtk file format. 

 ``radial_distribution`` : block 
    if set, the radial distribution function, g(r), is computed. The distribution represents the probability density
    of finding an atom at a given distance from a reference atom. It is computed by counting all atoms that are found 
    at a distance (r, r+dr) and normalized to the total volume and total number of atoms. Example

  radial_distribution
  {
     species = (In, Ga )
     cutoff_radius = 1.0 
     resolution = 0.001 
  }

  ``cutoff_radius`` : real
    is used to set the maximum distance of the distribution.

  ``resolution`` : real
    is used to specify the binning distance step (dr).

  ``species`` : specie names
    For each specie a radial distribution function is computed. The distribution is also resolved into the 
    other species present in the atomistic structure. Files are produced in output with names,
    structure_name_Specie_radial_suffix.dat, with column structure
    r   specie_1  specie_2   ...
   ...    ...       ...      ...

    In practice each column represents the probability density to finding a specie_i at distance r from the 
    reference specie. 
    The radial distribution can be computed also after VFF optimization by adding **radial_distribution** 
    to the plot line.





Importing an external file
^^^^^^^^^^^^^^^^^^^^^^^

An  external atomic  structure in a  valid  file  format may  be loaded  through the  keyword  *load_structure*.
Valid file formats are *xyz*, *gen* and  *tgn*. While *xyz* and  *gen* are commonly  used atomic descriptions, *tgn* format is an  internal format  used  in  tibercad which  contains information on  the  phisical regions  and mesh  elements associated to the atoms.

When  loading an external  structure file in  *xyz* format, the supercell vectors of  the  structure must  be  defined through the keyword *lattice_vectors*.
When   loading a gen  file this information is  already included.
After loading,  the external atomic  structure   is trimmed outside  the physical regions defined in  **Atomistic**  block. Then atoms  are  associated to  the current  mesh. 

**Note**: It  is  assumed that the atomic positions of  the  structure  are  consistent with the mesh.





Applying   periodicity
^^^^^^^^^^^^^^^^^^^^^^^

Periodicity may  be  applied  to  the  generated atomistic  structure in  1, 2  and  3D. 
Default  treatment of  periodicity is  different  for  each  mesh dimension,  to  reflect most  common cases.


For a  1D mesh, a  periodicity is  always applied in *y* and  *z*  directions. 
As  a  default, a  minimal cell  is  built in  *yz* plane (*x* being  the predefined growth direction). This   reflects the case  of  a  quantum well structure. The periodicity  vectors are  by  default those of the  minimal  cell  for  the  given material. It is  possible  to build   a  supercell in  the  yz plane,  by  defining the supercell  size through the **y_length** and **z_length**  keywords. In  this  case, appropriate  periodicity vectors will be  obtained and applied. It  is  possible to  add  a  periodicity along the *x* (growth) direction  with the  keyword  *periodicity = true*,   to design a superlattice structure.

For a  2D mesh, the  generated structure is  periodical  along  *z* axis. It  is  possible to  add  a  periodicity in *x* and *y* directions  with  the  keyword *periodicity = (<x-per>,<y-per>)*, where <x-per> and/or <y-per>  may  be  *true* or  *false*.

For a  3D mesh, the  generated structure is  not periodical  along  any axis (*cluster*).It  is  possible to  add  a  periodicity in any directions  with  the  keyword *periodicity = (<x-per>,<y-per>,<z-per> )*, where <x-per> and/or <y-per> and/or  <z-per> may  be  *true* or  *false*.

It must  be  noted that  here  we  refer  to  the  dimension of  the  geometrical  and  mesh models.
In  1 and  2D, however,  the atomic  stuctures are  in  any case built with at  least one layer along the  directions which are not  present in  the geometrical model, according to the material lattice parameters.




Use of  random alloy generator
^^^^^^^^^^^^^^^^^^^^^^^
By default, the atomic  structure of  an  alloy material is  generated based on   Virtual Crystal Approximation (**VCA**).
The VCA considers an alloy  :math:`{A_{x}B_{1-x}C}` as a fictitious material whose properties are a weighted average of the properties of the AC and BC components. In continuous media models, nonlinearities in the alloy parameters such as bandgap, effective masses or piezoelectric tensors are usually modelled with suitable bowing parameters, derived experimentally or from bulk DFT calculations. 
In  ETB,  tight-binding parameters are averaged between those of  the  alloy components,  weighting with the nominal alloy concentration.

An alternative  to VCA is the **random alloy** approach, where the  atoms of the components species are randomly distributed in the alloy region and their parameters are not  averaged.
Using atomistic approaches it is possible to model in detail realistic alloy  distributions within the nanostructures including the effect of local strain, thus leading to more accurate electronic calculations and allowing  to  study alloy fluctuactions in  a  region of a  device.


To  this  aim,  it  is necessary  to generate several random  samples of  the  alloy  structure,  in  order  to perform  a  statistical  analysis of  the  fluctuactions of  alloy  properties.
In TiberCAD, a  random  alloy structure is generated if  *random_alloy = true* in **Atomistic** block. In this case, the alloy  is created  by  random substitution of one of  the components of the host material with the substitution ion. For  example, InGaN alloy is  created by  substituting Ga atoms with In atoms in a GaN  lattice. The modality of  this  substitution is  ruled  by  *fix_mean_alloy_concentration*: if  this is  **true**, then a  fixed  number of  atoms (obtained  from  molar concentration)  is substituted, in randomly  chosen  positions of the  lattice. If *fix_mean_alloy_concentration* is **false**, then each Ga atom in the lattice is substituted with a probability proportional to the defined nominal molar fraction. This implies  that the  actual  number of  substitutions,  and  thus the  actual mean concentration of the  present  random  sample,  differs from  the nominal value. In  this  way  one  can  take  in  account  fluctuactions  of  the  global alloy  concentration.

Random generation is  ruled by  a  random seed,  which is  in  general different for  each  execution. In  case  one  needs  to  repeat calculations  with  exactly  the  same  random  configuration,  it  is  possible  to  keep it unmodified  by  defining the  same *random_generator_seed* in  all  the executions.


By default,  the  random alloy  structure is    built  with  an uniformly random approach. An option of  random  alloy  generation allows to  generate  **clustering** in  random  configuration. This  is  made  by  defining *clustering = true*. In  this  way,  the  alloy  structure is  not  built  entirely with  an uniformly random approach. Instead,  the  formation of clusters  of  substitution  atoms (e.g. In  atoms in  InGaN  alloy) is  favoured.  
The  keyword *clustering_seeds* define the fraction (default is 0.02) of  substitution ions which  are  placed  with an uniform distribution.  The  remaining fraction of ions  is placed with  a  higher  probability in  positions close to other  substitution ions.

Alloy  statistic, local  concentration...




Example
-----------------------


Here is an example of generation of  an  atomistic structure ::

  
 Region well1
  {
    material = InGaN
    x = 0.27
  }
 
  Region  qbarrierl
  {
    material = GaN
  }
   
  Region  qbarriert
  {
    material = GaN
 
  }
  

  

  Atomistic tb
  {
    reference_region = qbarrierl
    regions = (well1, qbarrierl, qbarriert)
    passivation = yes
   
    print = (xyz)

    clustering = false

    
    #load =  output_x0.20_1200/tb.xyz

    random_alloy = true
    random_generator_seed = 5
    plot_alloy_composition = true
    control_volume_radius = 0.5
    extract_alloy_statistics = true

    #z_length = 20 
    #y_length = 60  


  }




An  atomistic structure named *tb* is created based on  the  physical regions *well1, qbarrierl, qbarriert*.
The  reference region is *qbarrierl*,  which means that  the  reference lattice  will  be that  of  GaN.
Passivation of  surface  states  is  *true*. 
Since the  mesh  dimension is 1D, periodicity is  imposed along  y and  z  axis.

*random_alloy = true*  means  that  the  InGaN alloy in  the  well1  region  will  be generated  with  random  substitution.  The  random  distribution  will be  uniform,  since *clustering = false*.
We  define a *random_generator_seed = 5*,  which  can  be  used  in  a  following calculation  to  obtain  exactly  the  same  atomic  structure.
*extract_alloy_statistics = true* will  write info on  local concentration, obtained averaging on local  spheres defined by *control_volume_radius = 0.5*.
 
*plot_alloy_composition = true*  will plot the local concentration on  the  mesh  for  a  nice  visualization.







>>>>>>> .r4032
Module empirical_tb 
-----------------------


Module options
^^^^^^^^^^^^^^^^^^^^^^^


The following options influence the behaviour of the **Module**  ``empirical_tb``:

 ``name`` : string 
     defines the name of empirical_tb simulation

 ``regions`` : string 
    defines the regions associated to this ETB simulation

 ``atomistic_structure`` : string  
    defines the name of the atomistic structure which  will  be  used  by ETB  to  create the  TB Hamiltonian. It  must be  the name  of  one  structure created  through an **Atomistic** block in  **Device**  section. 
   

 ``potential_simulation`` : string
    defines the name of the simulation (e.g.  ``driftdiffusion``) that can provide electric potential

 ``strain_simulation`` : string 
    defines the name of the simulation (e.g. ``elasticity``) that can provide elastic strain



 ``Harrison_scaling`` : boolean
    if  **true** , scaling  of  ETB  parameters is  applied (usually in  presence  of  material  deformation which  causes atom  displacement  from  equilibrium  position).
    This  option is  by deafult  set to **true**   whenever  a  strain  simulation is  performed on  the  system, otherwise is  by  default  **false**.










Solver section
^^^^^^^^^^^^^^^^^^^^^^^


The Solver section of the **Module**  ``empirical_tb`` contains the following options:

 ``num_valence_eigenvalues`` : integer 
    defines the number of valence band eigenvalues and eigenfunctions to be found.

 ``num_conduction_eigenvalues`` : integer 
    defines the number of conduction band eigenvalues and eigenfunctions to be found.

 ``load_states`` : boolean 
    if **true**, a list of quantum states (eigenvalues  and eigenfunctions) are  loaded from  the  file specified by  ``load_path``. If the number of loaded electron or hole states is  lower  of  respectively *num_conduction_eigenvalues* or  *num_valence_eigenvalues*, then the remaining  states are calculated. 

 ``load_path`` : string 
    path of  the file  from  which  tb  states  can  be loaded

 ``long_tolerance`` : double 
    tolerance...

 ``guess_conduction`` : double
     the algorithm try  to  find  the closest eigenvalues to this absolute  value  of  energy (eV). If not defined, then by  default it is  calculated internally based on  the band edges.

 ``guess_valence`` : double
      the algorithm try  to  find  the closest eigenvalues to this absolute  value  of  energy (eV).If not defined, then by  default it is  calculated internally based on  the band edges. 


With  the  keyword *load_states = true* it is  possible to resume a  calculation starting from the last of the states loaded from a  file. Also,  one  can use  a  set  of already calculated states for an optical  spectrum calculation (see in the following). 
Just  run *empirical_tb simulation* to  load  the  states  and  then  *opticstb* to calculate spectrum  based on those   states.

  



Output
^^^^^^^^^^^^^^^^^^^^^^^

The available output variables for **Module**  ``empirical_tb`` , to be specified in the plot option, are the following:

 ``tbstates`` :
   states eigenfuctions in  *cube*  format

 ``MeshStatesNodes``: 
    plots  states  on  the  mesh. MeshStatesNodes is not supported in 1D calculations 

 ``MeshStates``: 
    plots  states  on  the  mesh in 1D calculations (ee.g. Quantum wells).  



By  default  the list of  eigenvalues[eV] is printed  in the output  directory in  the  output file *simulation_name*.dat.
This  file contains the  values of all  the  electron  and  hole  calculated eigenvalues (*EigenEnergy*), the Fermi level energy (*FermiLevel*) and  the  occupation index (*Occupation*) of  each  state. 

Through the keyword *tbstates*, the  states eigenfuctions may be stored  in  a  cube file for  visualization of wawefunctions with the open-source viewer **jmol** (*http://jmol.sourceforge.net*).

By defining *MeshStatesNodes* (or *MeshStates* in 1D), square  modules of  the  eigenstate wawefunctions  may be  plotted on  the FEM  mesh through the  output  file *simulation_name_msh*.vtu, e.g  tb_msh.vtu.



Module opticstb
----------------------

The **Module** ``opticstb`` implements the  calculation of optical properties based on  the  tb states; in general, the optical  matrix elements are calculated from the tb models specified in *initial_state_model* and *final_state_model*. In  the  current  version the  model must be  the  same  in  both cases. For  example ::

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
    *optical_spectrum_k_0* to  select  optical emission spectrum in *k=0* (:math:`\Gamma` point) 
    *matrix_elements* to  select calculation of  the optical matrix elements 


The output  file  for optical emission spectrum in *k=0* (:math:`\Gamma` point) is named *<simulation_name>_spectrum_k_0.dat*  and is composed by  several  columns  with the  following format ::


  energy [eV]  spontaneous_power_density_Px[W/eV] spontaneous_power_density_Py[W/eV]   
   spontaneous_power_density_Pz[W/eV] stimulated_power_density_Px[W/eV]     
  stimulated_power_density_Py[W/eV] stimulated_power_density_Pz[W/eV]    
   gain_Px[]  gain_Py[] gain_Pz[]


The output  file  for optical matrix elements is  named *<simulation_name>.dat* and  is  composed  by 5 columns with   the  following  format ::

  initial_state final_state |Px|^2 |Py|^2 |Pz|^2




Integrated spectrum
^^^^^^^^^^^^^^^^^^^^^^^


.....................




Example 1
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



 
Example 2
-----------------------

In  this  second  example,  we  will  see  how to  perform two  atomistic calculations on atomistic structure with different  size.
First  a  full  atomistic  structure will  be  generated,  representing  a quantum  dot with its  WL and buffer  layers. On  this full  structure,  a  VFF  relaxation will  be  executed. Then, EETB  will  be  applied  to  a  smaller subset  of  the  relaxed  structure, including only qdot  and WL.

First,  the  Device  structure  is  defined ::

  Device qdot
  { meshfile = qdot_Zb.msh
   mesh_units = 1e-9
   dimension = 3

   structure = zb
   x-growth-direction = (1,0,0)
   y-growth-direction = (0,1,0)
   z-growth-direction = (0,0,1)

   Region box
   {
     material = GaAs
   }

   Region buffer_up
   {
     material =  GaAs
   }

   Region buffer_down
   {
     material =  GaAs
   }

   Region Atomistic_cube_up
   {
     material =  GaAs
   }

   Region Atomistic_cube_down
   {
     material =  GaAs
   }

   Region dot
   {
     material = InGaAs 
     x = 0.7 
   }

   Region buffer_WL
   {
     material = InGaAs 
     x = 0.7 
   }

   Region Atomistic_cube_WL
   {
     material = InGaAs 
     x = 0.7 
   }


   # group together the regions for the quantum simulation
   Cluster Quantum_Atomistic
   {
     regions = (Atomistic_cube_WL,Atomistic_cube_down,Atomistic_cube_up,dot)
   }

   Cluster Quantum_buffer
   {
     regions = (Atomistic_cube_WL,Atomistic_cube_down,Atomistic_cube_up,dot,buffer_up,buffer_down,buffer_WL)
   }



Here,  two  **Clusters**  are  defined: *Quantum_buffer* refers to the full device  region where  we  want  to  generate an  atomistic description;  *Quantum_Atomistic* is  a  subset of  the previous one, and  will be  used  for  ETB calculations.



Next, the Atomistic block defines the atomistic structure *dot_atoms*, associated to **Cluster** *Quantum_buffer*  ::
   
   Atomistic dot_atoms
   {
   #regions specifica dove definire gli atomi
   regions = Quantum_buffer
   #reference_region stabilisce il reticolo di riferimento
   reference_region = Atomistic_cube_down
   #Enable hydrogen passivation
   passivation = yes
   #salva la struttura in vari formati, xyz per controllare con jmol, tgn per ricaricarla in tibercad
   #Nota: la struttura puo' essere ricaricata dal file tgn per risparmiare tempo
  #ma la MESH deve essere la stessa di quando la struttura e' stata salvata, altrimenti
   #i risultati non sono consistenti! Il file tgn viene salvato in path, l'xyz tra gli output
   print = (xyz, tgn)
   } 
  }


Next, Elasticity Module is defined. In  this  way,  macroscopic strain  is  calculated, deformation is  applied  to the  mesh (*mesh_deformation = true*) and displacement are  projected  to the  atomistic  structure (*strain_atomistic_structure = dot_atoms*).  ::


  Module elasticity
  {

  name = strain
  plot = (Strain,Stress,Displacement) 
  
  mesh_deformation = true
  shape_error = 0.1
  #shape_iterations = 3
  strain_atomistic_structure = dot_atoms
  
  Physics {
    body_force lattice_mismatch
    {
      reference_material = GaAs 
    }
  }

  Solver
  {
    relative_tolerance = 1e-6
  }

  Contact bottom {type = clamp}

  } 



The  result  of  **elasticity** is  taken as  a  initial  guess  for  next  calculations.  With  **vff** Module,atomic  positions  of  *dot_atoms* structure are  relaxed by  imposing  a  fixed  position for  all the  atoms  on  the  region  boundary (*boundary_conditions = all_around*)  ::


  Module vff
  {
   atomistic_structure = dot_atoms 
   boundary_conditions = all_around  
    plot = (xyz)
  }


Drift-diffusion Module will be  executed to  solve  Poisson at  equilibrium  ::


  Module driftdiffusion
  { 

  name = dd
  plot = (Ec, Ev, ElPotential, ElField, Polarization)
  regions = all
  coupling = poisson
  integration_order = 2

  #save_state = true 
  load_state = ./prova2/dd.tsv

  Physics
   {
    strain_simulation = strain
    recombination srh { }
    polarization (piezo, pyro) { } 
   }
  }


Finally,  we  define  **empirical_tb** Module  ::


  Module empirical_tb
  {
   
   name = tb
  
   regions = Quantum_Atomistic
 
   atomistic_structure = dot_atoms
  
   Harrison_scaling = false 
 
   potential_simulation = dd
  
   sparse_format = full
 
   #Plot in uscita, tbstates indica gli stati in formato cub o jvxl (plot sugli atomi)
   plot = (MeshStates,MeshStatesNodes)

  
   jmol_output_format = cube 

   Solver 
   {
    
    load_states = true
    load_path = prova

    num_valence_eigenvalues = 0 
    num_conduction_eigenvalues = 2 
    long_tolerance = 1e-7
   #PARAMETRI IMPORTANTI: il guess vicino al valore dell'autovalore fa si 
   #che non si prendano stati spuri dati dal folding dell'Hamiltoniana
   #internamente sono calcolati dei guess plausibili, ma per dot piccoli
   #lo shift e' tale che e' difficile dare un guess automatico. 
   #Il riferimento per il guess e' edge di valenza a 0.0
    guess_conduction = 0.15
    guess_valence = -0.10
   }
  }


Note  that ETB  will be  applied to the **atomistic_structure** *dot_atoms* (*atomistic_structure = dot_atoms*) , but only  to  the subset  of it whose  atoms  are  associated to the  regions included in **Cluster** *Quantum_Atomistic* (*regions = Quantum_Atomistic*).  ::



  Simulation
  {

   temperature = 300
   solve =  (strain,vff,dd,tb)
   resultpath = output 
   _format = vtk
  }


We  solve for  *strain* from **elasticity**,  then  for  *vff*  relaxation.  *dd*  calculates  band structure at equilibrium based on which *tb*  simulation is  applied.




