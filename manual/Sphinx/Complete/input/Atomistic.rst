

.. _Atomisticgen:


Atomistic  generator
=================================




The **Atomistic  generator** is  a   tool  of  tibercad which  allows  to  generate an atomic  basis associated  with the finite element mesh which belongs to  a  given  physical region, based  on  the material specifications and the growth directions defined  for  that region.
 
For  example, the  atomistic  structure ``tb``  is created  by  the  **Atomistic  generator** through the  following block of instructions contained in  the  *Device*  section  ::

  Atomistic tb
	{
	  reference_region = barrier_left
	  regions = Quantum 
	  passivation = yes
	  print = (xyz, gen, tgn)
          random_alloy = true

          radial_distribution
           { 
             species = (In,Ga)
             cutoff_radius = 1.0
             resolution = 0.001
           }

          alloy_statistics
           {
           control_volume_radius = 1.0
           plot_alloy_composition = true
           }

	}


The internal atomistic generator has the following features

* It supports any crystal structure with fcc, bcc, cubic and hexagonal Bravais lattice.
  Information about basis and primitive vectors is taken from the material database.

* It performs hydrogen passivation for any crystal structure supported. 

* For the  generation  of heterostructures  a pseudomorphic growth is assumed, with the same lattice constants of the  defined reference region.



The  definition of  a  *reference region* is  needed. It provides the  specifications on the material and parameters (growth direction, molar fraction, crystal phase) which  will be used to build the crystal structure. 

In this way,  only the atomic species are changed according to the different materials specifications, while the  lattice is  the reference material lattice. 
 
The structure may  be eventually relaxed by  projecting strain  calculated with a continuum elasticity model, or by  minimizing  the  energy  of  the  system applying  a  valence force field (VFF) model.


General Options
-----------------------



The  following kewords may  be  defined in **Atomistic**  block


 ``regions`` : string
    defines the Phisical Regions or  Cluster where to  build  the  atomistic  structure

 ``reference_region`` : string 
    defines the name of the region whose  material  and  growth directions  are  taken  as  a  reference for  lattice  construction

 ``reference_material`` : string 
    defines a general reference  material  and  growth directions that are taken as a  reference for  lattice  construction (may be useful if one chooses a  reference material not present in  any device region)

 ``translation`` : real vector
    defines a translation vector as *(dx, dy, dz)* used to move the atomistic structure w.r.t. the mesh regions. 
    It is useful for fine alignments of atoms to the mesh. 

 ``passivation`` : boolean
    if true an hydrogen passivation is  applied to the structure. 
    Hydrogen atoms are inserted at a defined distance in the same direction of the original crystal bond.

 ``hydrogen_distance`` : real
    when *passivation = true* , defines the bond length at which Hydrogen atoms  are placed from the surface atoms (default 1.2 Angstrom)

 ``periodicity`` : list of boolean
    add periodicity along one or more  of the three  directions *x,y,z*, depending on  the mesh dimension: for a 1D mesh, add periodicity in *x* direction (*periodicity = true*); for a 2D mesh, add periodicity along *x* and/or *y*  directions (e.g. *periodicity = (true,true)*; for a 3D mesh, add periodicity along *x,y* and *z* (e.g. *periodicity = (true,true,true)*. More examples in  the  following.

 ``supercell_size_y`` : real
    defines a lateral supercell size on *y-axis* (in Angstrom) in 1D (quantum wells) and 2D (nanowires) structures.

 ``supercell_size_z`` : real
    defines a lateral supercell size on *z-axis* (in Angstrom) in 1D and 2D structures.
 
 ``load_structure`` : string
    specifies a (relative) path/file for the loading of an external structure. Valid formats are *tgn, xyz, gen*.
    The external structure is trimmed outside  the regions on which the atomistic structure is defined in the input file. 
    Note: when loading *xyz* files, the supercell vectors must be defined  through the keyword *lattice_vectors* (see below)

 ``lattice_vectors`` : list of floating
    specifies the periodicity  vectors used when importing an external atomic structure in *xyz* format (see *load_structure*). The format is (x1,y1,z1,x2,y2,z2, x3,y3,z3),  values  are  in  A, e.g.( 52.728 , 0.0, 0.0 , 0.0, 45.664 , 0.0, 0,0, 42.758 ). 



 ``random_alloy`` : boolean
    if **true**,  a random alloy structure is generated (default = **false**).More examples in the following. 

 ``random_generator_seed`` : integer
    an integer defining the random seed. Useful to produce exactly the same random alloy structure in a  set  of  executions.   

 ``fix_mean_alloy_concentration`` : boolean
    if **true** the number of substituted atoms in a random alloy is as close as possible to the nominal molar fraction of the alloy.
    If **false** each atom is substituted with a probability proportional to the defined nominal molar fraction. The resulting actual number of substituted atoms may differ from the nominal  concentration value.

 ``clustering`` : specie name 
    Alloy stuctures  with clustering on the given atomic specie are produced, e.g. *clustering = In*. This is done by placing a fraction of seed 
    atoms (cluster_seed) with uniform distribution, then the remaing fraction of atoms is placed biased with a larger 
    probability near existing ions.  

 ``cluster_seeds`` : real
    defines the fraction of ions placed with uniform distribution. Must be > 0. (*default = 0.02*). 
     
 ``reference_atom`` : integer
    sets the atom label (typically 1=cations, 2=anions) where local concentration is computed. Default is cations.    

 ``meshdata_format`` : string
    defines the output format of the alloy statistics (default = "vtk") plotted when *plot_alloy_composition = true* (see next paragraph).

 ``print`` : string
    print the  generated atomic  structure with  one  of  the  following formats: *xyz*:  xyz file  format; *gen*: gen file  format; *tgn*: tibercad format which combines informations on the mesh and the atoms.



Optional Blocks
-----------------------



Two  optional blocks may  be  present in  **Atomistic**  block:


Block  **alloy_statistics**
^^^^^^^^^^^^^^^^^^^^^^^

 
This  block performs  the  calculation of   the local concentration  by averaging over small 
local spheres whose radius is defined by the  keyword *control_volume_radius*. 

Example:  ::

  alloy_statistics
  {
    control_volume_radius = 1.0 
    plot_alloy_composition = true
  }

The  following kewords may  be  defined in **alloy_statistics**  block


  ``control_volume_radius``: real
    sets the radius of the control volume ( *nm*), that is  of a sphere on which local concentration is computed (*default = 0.5 nm*).

  ``plot_alloy_composition`` : boolean
   if **true** the alloy composition is mapped on the mesh and plotted in a *vtk* file format. 


Block **radial_distribution**
^^^^^^^^^^^^^^^^^^^^^^^


If this block is declared, the radial distribution function, g(r), is computed. The distribution represents the probability density 
of finding an atom at a given distance from a reference atom. It is computed by counting all atoms that are found 
at a distance (r, r+dr) and normalized to the total volume and total number of atoms. 

Example:  ::

  radial_distribution
  {
     species = (In, Ga)
     cutoff_radius = 1.0 
     resolution = 0.001 
  }

The  following kewords may  be  defined in **radial_distribution**  block:


  ``cutoff_radius`` : real
    is used to set the maximum distance of the distribution.

  ``resolution`` : real
    is used to specify the binning distance step (dr).

  ``species`` : specie names
    For each specie a radial distribution function is computed. The distribution is also resolved into the 
    other species present in the atomistic structure. Files are produced in output with names given by
    *<structure_name_Specie_radial_suffix>*.dat, with column structure  ::

    r   specie_1  specie_2   ...
   

    In practice each column represents the probability density to find a specie *specie_i* 
    at a distance *r* from the reference specie. 
    The radial distribution can be computed also after VFF optimization, by adding **radial_distribution** 
    to the *plot* instruction.






Importing an external file
-----------------------



An  external atomic  structure in a  valid  file  format may  be loaded  through the  keyword  *load_structure*.
Valid file formats are *xyz*, *gen* and  *tgn*. While *xyz* and  *gen* are commonly  used atomic descriptions, *tgn* format is an  internal format  used  in  tibercad which  contains information on  the  phisical regions  and mesh  elements associated to the atoms.

When  loading an external  structure file in  *xyz* format, the supercell vectors of  the  structure must  be  defined through the keyword *lattice_vectors*.
When   loading a *gen*  file this information is  already included.
After loading,  the external atomic  structure   is trimmed outside  the physical regions defined in  **Atomistic**  block. Then atoms  are  associated to  the current  mesh. 

**Note**: It  is  assumed that the atomic positions of  the  structure  are  consistent with the mesh.





Applying   periodicity
-----------------------



Periodicity may  be  applied  to  the  generated atomistic  structure in  1, 2  and  3D. 
Default  treatment of  periodicity is  different  for  each  mesh dimension,  to  take into account  most  common cases.


For a  1D mesh, a  periodicity is  always applied in *y* and  *z*  directions. 
As  a  default, a  minimal cell  is  built in  *yz* plane (*x* being  the predefined growth direction). This   refers to  the case  of  a  quantum well structure. The periodicity  vectors are  by  default those of the  minimal  cell  for  the  given material. It is  possible  to build   a  supercell in  the  *yz* plane,  by  defining the supercell  size through the **supercell_size_y** and **supercell_size_z**  keywords. In  this  case, appropriate  periodicity vectors will be  obtained and applied. It  is  possible to  add  a  periodicity along the *x* (growth) direction  with the  keyword  *periodicity = true*,   to design a superlattice structure.

For a  2D mesh, the  generated structure is periodical  along  *z* axis. It  is  possible to  add  a  periodicity in *x* and *y* directions  with  the  keyword *periodicity = (<x-per>,<y-per>)*, where <x-per> and/or <y-per>  may  be  *true* or  *false*.

For a  3D mesh, the  generated structure is  not periodical  along  any axis (*cluster*). It  is  however possible to  add  a  periodicity in any directions  with  the  keyword *periodicity = (<x-per>,<y-per>,<z-per> )*, where <x-per> and/or <y-per> and/or  <z-per> may  be  *true* or  *false*.

It must  be  noted that  here  we  refer  to  the  dimension of  the  geometrical  and  mesh models.
In  1 and  2D, however,  the atomic  stuctures are  in  any case built with at  least one layer along the  directions which are not  present in  the geometrical model, according to the material lattice parameters.




Use of  random alloy generator
-----------------------



By default, the atomic  structure of  an  alloy material is  generated based on   Virtual Crystal Approximation (**VCA**).
The VCA considers an alloy  :math:`{A_{x}B_{1-x}C}` as a fictitious material whose properties are a weighted average of the properties of the AC and BC components. In continuous media models, nonlinearities in the alloy parameters such as bandgap, effective masses or piezoelectric tensors are usually modelled with suitable bowing parameters, derived experimentally or from bulk DFT calculations. 
In  ETB,  tight-binding parameters are averaged between those of  the  alloy components,  weighting with the nominal alloy concentration.

An alternative  to VCA is the **random alloy** approach, where the  atoms of the components species are randomly distributed in the alloy region and their parameters are not  averaged.
Using atomistic approaches it is possible to model in detail realistic alloy  distributions within the nanostructures including the effect of local strain, thus leading to more accurate electronic calculations and allowing  to  study alloy fluctuactions in  a  region of a  device.


To  this  aim,  it  is necessary  to generate several random  samples of  the  alloy  structure,  in  order  to perform  a  statistical  analysis of  the  fluctuactions of  alloy  properties.
In TiberCAD, a  random  alloy structure is generated if  *random_alloy = true* in **Atomistic** block. In this case, the alloy  is created  by  random substitution of one of  the components of the host material with the substitution ion. For  example, InGaN alloy is  created by  substituting Ga atoms with In atoms in a GaN  lattice. The modality of  this  substitution is  ruled  by  *fix_mean_alloy_concentration*: if  this is  **true**, then a  fixed  number of  atoms (obtained  from  molar concentration)  is substituted, in randomly  chosen  positions of the  lattice. If *fix_mean_alloy_concentration* is **false**, then each Ga atom in the lattice is substituted with a probability proportional to the defined nominal molar fraction. This implies  that the  actual  number of  substitutions,  and  thus the  actual mean concentration of the  present  random  sample,  differs from  the nominal value. In  this  way  one  can  take  into  account  fluctuactions  of  the  global alloy  concentration.

Random generation is  ruled by  a  random seed,  which is  in  general different for  each  execution. In  case  one  needs  to  repeat calculations  with  exactly  the  same  random  configuration,  it  is  possible  to  keep it unmodified  by  defining the  same *random_generator_seed* in  all  the executions.


By default,  the  random alloy  structure is    built  with  an uniformly random approach. An option of  random  alloy  generation allows to  generate  **clustering** in  random  configuration. This  is  made  by  defining the atomic specie subjected to clustering with *clustering = specie_name*, e.g. *clustering = In* (by default no clustering is applied). In  this  way,  the  alloy  structure is  not  built  entirely with  an uniformly random approach. Instead,  the  formation of clusters  of  substitution  atoms (e.g. In  atoms in  InGaN  alloy) is  favoured.  
The  keyword *clustering_seeds* define the fraction (default is 0.02) of  substitution ions which  are  placed  with an uniform distribution.  The  remaining fraction of ions  is placed with  a  higher  probability in  positions close to other  substitution ions.






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

 #   clustering = In

    
    #load =  output_x0.20_1200/tb.xyz

    random_alloy = true
    random_generator_seed = 5
    plot_alloy_composition = true
    control_volume_radius = 0.5
    extract_alloy_statistics = true

    #supercell_size_z = 20 
    #supercell_size_y = 60  


  }




An  atomistic structure named *tb* is created based on  the  physical regions *well1, qbarrierl, qbarriert*.
The  reference region is *qbarrierl*,  which means that  the  reference lattice  will  be that  of  GaN.
Passivation of  surface  states  is  *true*. 
Since the  mesh  dimension is 1D, periodicity is  imposed along  y and  z  axis.

*random_alloy = true*  means  that  the  InGaN alloy in  the  well1  region  will  be generated  with  random  substitution.  The  random  distribution  will be  uniform,  since *clustering = In* is commented.
We  define a *random_generator_seed = 5*,  which  can  be  used  in  a  following calculation  to  obtain  exactly  the  same  atomic  structure.
*extract_alloy_statistics = true* will  write info on  local concentration, obtained averaging on local  spheres defined by *control_volume_radius = 0.5*.
 
*plot_alloy_composition = true*  will plot the local concentration on  the  mesh  for  a  nice  visualization.



