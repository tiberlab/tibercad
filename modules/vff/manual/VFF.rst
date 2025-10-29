..   <marker>

.. _VFFTheory:


Valence Force Field
=================================

Theory
--------------------------------

The Valence Force Field (VFF) is an empirical method which allows to calculate the total energy of crystaline structures as a sum of two center and three center terms which depend on bond stretchind and bond bending. The total potential :math:`U` is calculated as a sum of atomic contributions :math:`U=\sum_{i}U_i` where :math:`i` is an index running on the atoms and the form of :math:`U_{i}` depends on the crystal structure.TiberCAD includes potentials suitable for Zinc Blende and Wurtzite crystals. 

In the case of Zinc Blende the original Keating potential is used:

.. math::

    U_i =\frac{3\alpha}{16 r_{0}^2}\sum_{j=1}^{4} \left(\textbf{r}_{ij}-r_{0}^{2} \right)^{2} + \frac{3\beta}{8 r_{0}^2}\sum_{j=1}^{4}\sum_{k>j}^{4} \left(\textbf{r}_{ij}\textbf{r}_{ik}-r_{0}^{2}\cos \theta_{0} \right)^{2} 

where :math:`\alpha` and :math:`\beta` are the bond stretching and bond bending constants,  :math:`r_{0}`  and :math:`\theta_{0}` the reference bond length and bond angle of the bulk crystal. 

In the case of Wurtzite a modified Keating potential is used. In Zinc Blende each atom is bonded by four identical bonds to four first neighbor atoms placed on the vertices of an ideal tetrahedorn. In Wurtzite crystals the neirest neighbors are placed on the vertices of a distorted tetrahedron, such that three bonds are identical but the fourth one has a different length an angle. The non ideality is usually described by an internal parameter :math:`u`. 

Therefore in the case of Wurtzite the following potential is used:

.. math::

    U_{i} = & \frac{3\alpha}{16 r_{0}^2}\sum_{j=1}^{3} \left(\textbf{r}_{ij}-r_{0}^{2} \right)^{2} + \frac{3\alpha '}{16 r_{0}'^2} \left(\textbf{r}_{i4}-r_{0}'^{2} \right)^{2} + \\
            &+ \sum_{j=1}^{3}\sum_{k>j}^{3}\frac{3\beta}{8r_{0}^2} \left(\textbf{r}_{ij}\textbf{r}_{ik}-r_{0}^{2}\cos \theta_{0} \right)^{2}+\frac{3\beta '}{8r_{0}r_{0}'}\sum_{k=1}^{3} \left(\textbf{r}_{i4}\textbf{r}_{ik}-r_{0}r_{0}'\cos \theta_{0}' \right)^{2}

where different reference bond length and angles are used if bonds along the :math:`[0001]` crystal direction are considered. In principle, also differen coefficients :math:`\alpha`, :math:`\alpha '` and  :math:`\beta`, :math:`\beta '` can be used. Once the potential as a function of atomic coordinates is known, the system geometry can be relaxed by minimizing the VFF potential. 

The stiffness constants used in continuum elasticity can be written in terms of VFF coefficients as well, for further details refer to the work of Keating [P.N. Keating Phys. Rev., 145 (1966), p. 637] for the Zinc Blende and to the work of Camacho and Niquet [D. Camacho and Y.M.Niquet, Physica E, 42 (2009), p.1361], even though through a set of linearly dependent equations. Therefore it is not possible to defin coefficients which correspond to a given set of stiffness constants. If no parameters are provided via input file by the user, tiberCAD will calculate the parameters internally by minimizing the mismatch between the stiffness constants provided in database and the one calculated from the VFF. This allows to ensure the best consistency possible between continuum and atomistic model.





Module vff
-----------------------


In  tiberCAD,  it  is  possible  to  perform a relaxation of  an atomistic structure by  means of 
the **Module**  ``vff``. Here is  an example  :: 


  Module vff
  {
   regions = all
   plot = (StrainNodes) 

   atomistic_structure = tb
   boundary_conditions = all_around
 
   Solver
   {
    absolute_tolerance = 1e-4
    print_level = 0 
   } 
  }


A VFF calculation is desirable in combination with other atomistic methods, such as empirical tight binding, as it allows to take into account internal strain.



Module options
^^^^^^^^^^^^^^^^^^^^^^^



The VFF module is controlled by the following options:

``atomistic_structure`` : string
  Defines the atomistic structure on which the VFF must be calculated. The atomistic structure must have been declared in the appropriare Device section.

``boundary_conditions`` : string
  Defines if and how some atoms are fixed to simulate different boundary conditions. The available options are ``substrate``, ``free_standing`` and ``all_around``. The default is ``free_standing``, meaning that all the atoms are allowed to relax. If ``all_around`` is specified, all the outer atoms are fixed. Each atom which is not bonded to 4 atoms which belong to the same structure, or each atom bonded to a passivation hydrogen, is considered an outer atom. If ``substrate`` si specified, all the atoms belonging to an initial layer orthogonal to a given direction are fixed. This last keywordis used to simulate epitaxial heterostructures. The behaviour of ``substrate`` is determined by additional keywords here enlisted.

``substrate_plane`` : string
  It can take values ``x``, ``y`` or ``z``. If ``boundary_conditions`` is set to ``substrate``, it specifies the direction orthogonal to the substrate in the cartesian reference. The atoms with smallest coordinates, within a given tolerance, along the given direction are fixed (i.e., it is assumed that the structure ir oriented along positive cartesian direction, with substrate underneath). If  ``boundary_conditions`` is not set to ``substrate``, this keyword has no effect.

``substrate_updown`` : bool
  If set to true, both atoms with smallest and largest coordinates are fixed, to simulate a subatret condition on both the lowest and highest structure plane. This can be useful to combine different relaxation methods. If  ``boundary_conditions`` is not set to ``substrate``, this keyword has no effect.

``substrate_tol`` : double
 Define the thickness of the layer of fixed atoms (in Angstrom). The default is ``1.0``, which is usually enough to fix the first atomic layer. If  ``boundary_conditions`` is not set to ``substrate``, this keyword has no effect.


Solver options
^^^^^^^^^^^^^^^^^^^^^^^


The VFF supports a structure optimization through a non-linear conjugate gradient algorithm. 
he tolerance is expressend as a force in ``eV/A`` and the default value is ``1e-3``. The tolerance can be modified with the keyword ``absolute_tolerance``. Example:

::

  Solver
  {
  absolute_tolerance = 1e-4
  }







Physics
^^^^^^^^^^^^^^^^^^^^^^^




The Physics section of the input file is used to allow the definition of a **custom** bond stretching parameter :math:`\alpha` and a **custom** bond bending parameter :math:`\beta` through special blocks **keating** with the type *user*. As a default, all these parameters are calculated internally. An example of custom parameters is ::

  Physics
  {
   keating user
   {
    regions = InN
     alpha = 67.38
     beta = 10.01
   }
   keating user
   {
    regions = GaN
    alpha = 88.35
    beta = 20.92
   }
  }

In this example the user explicitely set the model coefficients for all the regions made of Indium Nitride and for all the region made of Gallium Nitride. The coefficients are specified in ``N/m``. 

If at least one of the material is an **alloy**, the user may choose between a **random alloy** approximation or a **virtual crystal approximation (VCA)**. In a **random alloy** approximation the two center parameters are assigned according to a reference bulk material which depend on the bond atoms (e.g. a Ga-N bond will get bulk GaN parameters assigned), and three center parameters will be assigned as an average of corresponding bulk terms (e.g. a In-N-Ga bond will have parameters averaged between bulk InN and GaN). In the case of **random alloy** the user can specify both components as in the following example, where we consider two regions, one associated to alloy InGaN and one associated to GaN.   ::

  Physics
  {
   keating user
   {
    regions = InGaN
    component InN
    {
     alpha = 67.38
     beta = 10.01
    }
    component GaN
    {
     alpha = 88.35
     beta = 20.92
    }
   }
   keating user
   {
    regions = GaN
    alpha = 88.35
    beta = 20.92
   }
  }


For InGaN, parameters for each alloy component, InN and GaN are defined separately through the subblock **component**, where the keyword *component* is followed by the name of the alloy component, e.g. *InN*, as in the following ::

  component InN
   {
     alpha = 67.38
     beta = 10.01
   }



In the **virtual crystal approximation** the coefficients are automatically calculated in order to approximate properly the alloy stiffness coefficient calculated with Vegard law, in order to ensure consistency with continuum models. A user may override these value by specifying the alloy coefficients, as in this example :: 

  Physics
  {
   keating user
   {
     regions = InGaN
     alpha = 70.05
     beta = 12.00
   }
   keating user
   {
    regions = GaN
    alpha = 88.35
    beta = 20.92
   }
  }

Different coefficients can be defined for different regions just by specifying  the correct region names. Note that this settings have to be consistent with the definition ``random_alloy`` in the **Atomistic Generator**.   



Output
^^^^^^^^^^^^^^^^^^^^^^^




The available output variables for **Module**  ``vff`` , to be specified in the plot option, are the following and correspond to the  strain tensor extrapolated from aton positions in the unit cell after relaxation:

 ``StrainNodes`` :
   strain tensor associated to nodes

 ``StrainCells``: 
    strain tensor associated to cells  




Example
--------------------------

First,  the  Device  structure  is  defined for a InGaN/GaN quantum well heterostructure  ::

  Device 
  {
   meshfile = well.msh 
   mesh_units = 1e-9    #nanometers
   dimension = 1 

   structure = wz
   x-growth-direction = (0,0,0,1)
   y-growth-direction = (1,0,-1,0)
   z-growth-direction = (-1,2,-1,0)
 

   Region QW
   {
    material = InGaN
    x=0.1
   }

   Region barrier_left
   {
    material = GaN
   }

   Region barrier_right
   {
    material = GaN
   }

   Region device_left
   {
    material = GaN
   }

   Region device_right
   {
    material = GaN
   }

   # group together the regions for the quantum simulation

   Cluster Quantum
   {
    regions = (barrier_left, barrier_right, QW)
   }

   Atomistic tb
   {
    reference_region = device_left
    regions = Quantum 
    passivation = no 
    print = (xyb, xyz, gen, tgn)
   
   }

   }#Device

An  atomistic block defines the region where  an atomistic structure needs to be  generated, 
that is the  QW  and the barrier  regions.
Based on  the  GaN/InGan/GaN  atomistic  structure generated,  ETB  calculations will  be performed, so that we  define a **Module** ``empirical_tb``  ::


  Module empirical_tb
  {
  regions = Quantum
  name = tb
  atomistic_structure = tb
  Harrison_scaling = true 

  strain_simulation = strain
  potential_simulation =  driftdiffusion

  sparse_format = full   

  projection_length = 1.0
  dangling_bond_scaling = 1.0
  dangling_bond_onsite = -5.0

  plot = (tbstates, ProbabilityDensity)
  jmol_output_format = cube

  #assemble_hamiltonian = false 

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


A  **Module** ``opticstb`` is  defined to  calculate optical properties with  ETB  ::

  Module opticstb
  {

  name = spectrum 
  plot = optical_spectrum_k_0 

  initial_state_model = tb
  final_state_model = tb


  Emin = 3.0
  Emax  = 5.0
  dE = 0.001

  }


A  **Module**  ``elasticity`` is  defined to calculate macroscopic strain  as  a first  guess 
for  VFF relaxation  ::


  Module elasticity
  {

  name = strain
  mesh_deformation = true
 

  plot = (Strain,Stress,Displacement)

  Physics
  {

   body_force lattice_mismatch
   {
    
   reference_material = GaN
   x-growth-direction = (0,0,0,1)
   y-growth-direction = (1,0,-1,0)
   z-growth-direction = (-1,2,-1,0)

  }

  }

   Contact cathode{type = clamp}

  }

The  **Module** ``vff`` is  here defined in  such  a  way to be applied to  the  whole atomistic structure *tb*,  with boundary  conditions fixing all the outer atoms (*all_around*). 
In output, strain tensor extrapolated from vff  relaxation is plotted (*StrainNodes*)  ::


  Module vff
  {
   regions = all
   plot = (StrainNodes) 

   atomistic_structure = tb
   boundary_conditions = all_around
 
   Solver
   {
    absolute_tolerance = 1e-4
    print_level = 0 
   } 
  }

Finally, drift-diffusion is  defined for an  equilibrium calculation of the band profiles  ::

  Module driftdiffusion
  { 

  name = driftdiffusion
  regions = all
  coupling = poisson

  #plot = (Ec, Ev, eQFermi, hQFermi, ElField , eDensity, hDensity, Polarization)
  plot = (Polarization, Ec, Ev, eQFermi, hQFermi, ElField, eDensity, hDensity, ElPotential)

  save_state = true
  #load_state = dd.tsv
  # solve_after_load = true

  Physics
  {

     strain_simulation = strain

     use_density_predictor = true
     background_conductivity = 1e-6


     recombination   srh {
	}

     recombination  direct
     { 
        C = 1.1e-8
     }

    mobility
    { 
      type = doping_dependent
      #model = field_dependent
      #low_field_model = doping_dependent
    }
  } #endphysics 


  Contact cathode
  {
    type = ohmic
    voltage = $Vb[0.0]  
    rec_velocity_p = 1e7
  }

  Contact anode
  {
    type = ohmic
    voltage = 0.0
   rec_velocity_n = 1e7
  }

  } #driftdiffusion



We  solve first for **strain**  with  *elasticity*, which applies in this  case 
a  mesh  deformation to the  device, accounting for the  calculated displacements.
Then, based  on this deformed mesh,  **vff** is  applied,  which  relaxes the  atomistic  structure associated  to the  mesh.  Finally, an equilibrium  solution is  found for potentials and  then empirical tight-binding Hamiltonian is calculated on the  relaxed atomic  structure (**tb**), to  obtain  eigenstates and  optical properties (**spectrum**)  ::

 

  Simulation
  {
   verbose = 4
   temperature = 300

   solve = (strain, vff , driftdiffusion,tb,spectrum)

   resultpath =  output
   output_format = grace

  }


.. note:: 
          Note that it is always highly advisable to solve first elasticity, before a VFF relaxation. This will provide a reasonable first guess for the VFF solver and also will apply the correct deformation to the mesh, assuring that atoms are still associated to the correct FEM elements.

