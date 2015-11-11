..   <marker>

.. _ETBTheory:


Empirical Tight Binding
=================================


Theory
--------------------------------


Tight Binding approach
^^^^^^^^^^^^^^^^^^^^^^^

The Tight binding method (TB) is a technique to describe
the electronic and optical properties of a material through
an atomistic description. A corpuscolar approach is needed
in devices with nanometric features, where a detailed description
of electronic ad optical properties cannot be reached with an
effective description.

In the following we  will give a brief description of general features
of Tight Binding Model. Then we will focus on the  technique implemented
in TiberCAD, namely Empirical Tight Binding (ETB).


The TB model has been introduces by Slater and Koster in 
a notable publication in 1954
(see [SlaterKoster]_).

The starting point is the *LCAO*  (linear combination
of atomic orbitals) representation. 

If we assume that the electronic
states :math:`\Psi`
in the system under observation are perturbations
of the atomic states, it's convenient to use the LCAO basis:

.. math::
    :label: lcao_basis 

    \left|\Psi\right\rangle=\sum_{\alpha \Vec{R}}C_{\alpha\Vec{R}}\left|{\alpha,\Vec{R}}\right\rangle


Where *R* is the atom position and :math:`\alpha` a quantum number
to distinguish between different orbitals.



Using the previous relation  in the time independent Schroedinger equation,
we obtain  the expression:

.. math::
    :label: eq_schroed

    \sum_{\alpha'\Vec{R}^{'}}C_{\alpha'\Vec{R}^{'}}\left[H_{\alpha'\Vec{R}^{'}\alpha\Vec{R}}-ES_{\alpha'\Vec{R}^{'}\alpha\Vec{R}}\right] = 0


where 

:math:`H_{\alpha'\Vec{R}^{'}\alpha\Vec{R}}` and :math:`S_{\alpha'\Vec{R}^{'}\alpha\Vec{R}}`
are the hamiltonian matrix element and the overlap matrix and 
are respectively given by:

.. math::
    :label: matrix_elem

    H_{n'\alpha',n\alpha} = \left\langle n'\alpha'\right| H \left|n\alpha\right\rangle 


    S_{n'\alpha',n\alpha} = \left\langle n'\alpha'|n\alpha\right\rangle



The Hamiltonian may possibly  include an external potential 
:math:`V_{ext}` term, including an externally applied potential and mean field
corrections. 
The overlap matrix is the identity when the basis is orthogonal.
This is not the case in general, since usually  atomic orbitals have non zero
overlap. However,  the basis set can be orthogonalized using a 
Lowdin orthogonalization. Anyway, resulting wavefunctions
have longer range, and it determines  a worst numerical formulation
than when using a non orthogonal basis. 
The overlap matrix *S* can be calculated from the basis function, but 
in order to solve the generalized eigenvalue problem 
the matrix elements of the Hamiltonian
need to be evaluated. The method used to obtain these matrix elements
classifies the different TB implementations.


The Hamiltonian matrix elements can be classified in four categories
(see [DiCarloreview]_)


  * ``On-site``:  when the atomic wavefunctions and the potential are centerd on the same site

  * ``Two-centre``:  when wavefunction and potential are on the same site and the other wavefunctions are on different sites

  * ``Three-centre``:  when the wavefunctions and the potential are all on different sites

  * ``Local environment correction``:  when the two wavefunctions are on the same site and the potential is on a different site

 
In the following we will consider only on-site and two-centre 
contributions,  a reasonable trade-off between
computational effort and accuracy.


Empirical Tight Binding
^^^^^^^^^^^^^^^^^^^^^^^

As seen in previous section, the type of matrix element included 
in a TB representation is one of the discriminants between one model and
another.
The other features which differentiate  TB models
are the number of interacting neighbours (usually nearest neighbours or
second nearest neighbours), the orthogonality of the basis set and the method 
used to calculate the matrix elements.
Focusing on the last point, there are three main techniques to evaluate 
these parameters: empirical, semi-empirical and ab initio methods.

In empirical tight binding (ETB), matrix elements are calculated
as fitting parameters of characteristic bulk quantities, such as
effective mass, energy gap and split-off energy.


In order to
solve the Schrodinger equation 
we need to evaluate
two different kinds of matrix element: *on-site*  matrix elements and
*hopping*  matrix elements.

The **empirical tight binding** technique consists in finding these quantities
by fitting the bulk properties of semiconductors, assuming an orthogonal basis,
i.e. *S =* **I** with **I** the identity matrix. 


The physical meaning of the on-site term is straight: it's the eigenenergy of
the corresponding orbital. The hopping energy can be written in different notations.
In the original work of Slater and Koster (see [SlaterKoster]_), 
hopping elements are given by
linear combinations of atomic orbitals in a two centre approximation.
Note that Slater and Koster work is based on a Lowdin orthogonalized basis,
leading to zero overlap between orbitals. This is a crucial point in ETB
technique, as it allows  to solve an eigenvalue problem instead of a generalized
eigenvalue problem, with a much lower computational effort, and it's the key
for the efficiency of ETB respect to ab initio techniques.

It means that any interatomic parameter is labeled as :math:`V_{\alpha\beta\mu}`,
where :math:`\alpha, \beta = (s, p, d)` are the atomic orbitals involved and
:math:`\mu` is the component of angular momentum around the bond axis, 
i.e. :math:`\mu = \sigma`, :math:`\mu = (\sigma, \pi)`, :math:`\mu=(\sigma, \pi, \delta)`
respectively for :math:`ss, sp, sd, pp` and :math:`dd` bonds. 
This is the so called **molecular orbital** notation. 
The relationship between these quantities and hopping matrix elements
(Koster notation) is shown in Slater and Koster work.  



The chemestry of localized states in covalent semiconductors needs at 
least an eight band parameterization (one *s* and three *p* orbitals with spin
degeneracy) 
(see [Vogl]_ ).



Historically, the first approach to empirical
tight binding has been performed by using these :math:`sp^{3}` basis, often expanding the interaction to second nearest neighbours.
In fact first-neighbour :math:`sp^{3}` models fail in 
many points: it's proved that they can not
reproduce an indirect gap in diamond and zincblende materials
and they cannot fit even the lowest conduction bands
of semiconductors as Ge, Si, AlAs and GaP.

However a second-neighbour basis is not desirable as it
cannot be used to describe random alloys and heterointerfaces in a straight way,
while first neighbour parameterizations are used under the assumption
that parameters are portable when we're not dealing with bulk structures,
and locally the neighbours interact as in the corresponding bulk material.

It is possible to overcome these deficiencies by including 
an excited *s* state, i.e. the :math:`s^{*}` state, on each atom (see [Vogl]_). 
With this approach it's possible to get a ten band nearest neghbour
:math:`sp^{3}s^{*}` parameterization which describe with a good degree of 
accuracy bulk IV group semiconductors, zincblend and wurtzite III-V
materials and ternary alloys. 

The :math:`sp^{3}s^{*}`  parameterization has been used for years and it's still
routinely applied. However, it still suffers some inner 
limitations. For example, it cannot fit the *X* valley transverse mass
, making its application critical
in cases where *X* valley contribution are important, 
such as higly confined structures or materials with a band minimum in 
*X* like Silicon. There are two ways to overcome
these limitation: increasing the number of neighbours or
increasing the number of parameters.


For the reasons previously explained, we find
more convenient not to discard  the nearest neighbour approach.
Luckily we can rely on modern parameterizations
that include *d* orbitals, leading to a *30* bands description.
The :math:`sp^{3}s^{*}d` parameterization by Jancu 
(see [Jancu]_ , [Jancuwz]_ ), 
gives a very accurate
description of  C, Si, Ge, AlP, InP, GaAs, AlAs, InAs, GaSb, AlSb, InSb, 
GaN, AlN and InN. This parameterization introduces
seven additional hopping matrix elements but, at the
cost of an increased computational cost, it offers a very accurate description
of nanostructures. 







Module empirical_tb 
-----------------------


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




In the prevoius  example, TB  calculations  are applied  to  the  physical  region ``Quantum`` (usuallly a  Cluster) ::

  regions = Quantum   


An  atomistic structure must be  coupled with  this physical region, in  order  to  apply the TB  model ::

  atomistic_structure = tb

This  is  made through the **Atomistic  generator**,   the   tool  of  tibercad which  allows  to  generate an atomic  basis associated  with the finite element mesh which belongs to  a  given  physical region, based  on  the material specifications and the growth directions defined  for  that region.



See   section :ref:`Atomisticgen`  for  a  detailed description of the  Atomistic  generator.








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
    passivation = yes
    print = (gen, xyz, tgn)

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
    plot = (tbstates, MeshStatesNodes)

    Solver
     {
      num_valence_eigenvalues = 1 
      num_conduction_eigenvalues = 1
      long_tolerance = 1e-4  
      }

  }


Note  that  the  results of potential (from driftdiffusion)  and  strain (from  elasticity) simulations   are  recalled,  to  project  the  correct  potential  profile  and  atom  displacement  to  the  TB  Hamiltonian.
In  output,  the ground  states  of  conduction  and  valence  bands  are saved  together  with  the  wavefunctions.States  are  also  plotted  on  the  FEM  grid (through *MeshStatesNodes*).

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


The  output of  calculation is  composed  by  the *matrix elements* file *opt.dat* and  by the  optical spectrum in :math:`\Gamma` point, in  file *opt_spectrum_k_0.dat*.
For  the  latter output, the  keyword ::

  compute_strengths = true

is  necessary,  to calculate   emission spectrum power  assuming  all  occupations equal to  one,  since  the  simulation is  performed  at  equilibrium and  thus with nominal zero  occupation of states.
  




 
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
     regions = (Atomistic_cube_WL,Atomistic_cube_down,Atomistic_cube_up,dot,
       buffer_up,buffer_down,buffer_WL)
   }



Here,  two  **Clusters**  are  defined: *Quantum_buffer* refers to the full device  region where  we  want  to  generate an  atomistic description;  *Quantum_Atomistic* is  a  subset of  the previous one, and  will be  used  for  ETB calculations.



Next, the Atomistic block defines the atomistic structure *dot_atoms*, associated to **Cluster** *Quantum_buffer*  ::
   
   Atomistic dot_atoms
   {
   
   regions = Quantum_buffer
   #reference_region defines the reference lattice  
   reference_region = Atomistic_cube_down
   #Enable hydrogen passivation
   passivation = yes
  
   print = (xyz, tgn)
   # save the structure in  several formats 
   } 
  }


Next, Elasticity Module is defined. In  this  way,  macroscopic strain  is  calculated, deformation is  applied  to the  mesh (*mesh_deformation = true*) and displacement are  projected  to the  atomistic  structure (*strain_atomistic_structure = dot_atoms*).  ::


  Module elasticity
  {

  name = strain
  plot = (Strain,Stress,Displacement) 
  
  mesh_deformation = true
  shape_error = 0.1
 
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
  # load_state = ./prova2/dd.tsv

  Physics
   {
    band_properties 
     {
      density_of_states bulk_kp 
      {
        strain_simulation = strain
      }
     }
    
    recombination srh { }
    polarization (pyro) {}
    polarization (piezo) {strain_simulation = strain}
 
   }
  }

In  alternative,  a  previuosly calculated solution  from  drift-diffusion may be  loaded through ::

    load_state = <path>

Finally,  we  define  **empirical_tb** Module  ::


  Module empirical_tb
  {
   
   name = tb
  
   # Physical  regions to which etb is applied
   regions = Quantum_Atomistic
 
   # atomistic  structure on  which etb is calculated
   atomistic_structure = dot_atoms
  
   potential_simulation = dd
  
   plot = (MeshStates,MeshStatesNodes)

   jmol_output_format = cube 

   Solver 
   {
    
    #  load_states = true
    #  load_path = prova

    num_valence_eigenvalues = 0 
    num_conduction_eigenvalues = 2 
    long_tolerance = 1e-7

    guess_conduction = 0.15
    guess_valence = -0.10
   }
  }


By  default, reasonable  guesses are  generated  internally to  help  the solver  to get convergence. 
Usually these  values  are  good  enough to  get  a steady covergence  of  solutions.
However, in  this particular case of small  energy  gap  material (InGaAs), 
it  may be  necessary to  set manually a  guess close to the expected value of  the eigenstate. ::

   guess_conduction = 0.15
   guess_valence = -0.10

In this way, we  avoid erroneous states due to the  folding of the Hamiltonian.
The reference of the guess is the valence top band edge at 0.0 eV.


Note  that ETB  will be  applied to the **atomistic_structure** *dot_atoms* (*atomistic_structure = dot_atoms*) , but only  to  the subset  of it whose  atoms  are  associated to the  regions included in **Cluster** *Quantum_Atomistic* (*regions = Quantum_Atomistic*).  ::



  Simulation
  {

   temperature = 300
   solve =  (strain,vff,dd,tb)
   resultpath = output 
   format = vtk
  }


We  solve for  *strain* from **elasticity**,  then  for  *vff*  relaxation.  *dd*  calculates  band structure at equilibrium based on which *tb*  simulation is  applied.




