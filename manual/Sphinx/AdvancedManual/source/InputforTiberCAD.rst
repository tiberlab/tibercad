.. _InputforTiberCAD:

Input for TiberCAD
==================

Input for TIBERCAD is composed by an input file e.g. ``"input.tib"`` and a mesh file generated 
by a mesher software: as for now, mesh files from GMSH (*.msh, v.1 and v.2.0 )
and from ISE-TCAD (*.grd) are supported.

Be sure that the material files are in the correct directory ( as specified in  :ref:`SimulationSection` ).

To run the program, type: **tibercad** *input_file_name*


Description of Input file structure
---------------------------------------

A valid input file for TIBERCAD is a text file with the structure described in the following.
In the whole input file, everything following a '#' is considered as a comment and is
disregarded; blank lines can be present anywhere and are disregarded too.

Input file is composed by several **sections**: each **section** begins with a **sectionname** preceded by **"$"** (e.g. $Physics).

A section is enclosed between **"{"** and **"}"** brackets and is possibly composed (de-
pending on the section) by a variable number of **blocks** enclosed between "{" and "}"
brackets.

Each **block** can be possibly composed by one or more blocks, each preceded by a
block-name.

The elementary block (**parameters-block**) is a block which contain zero or any
number of **parameter assignements** in the form::

  "*tagname = tagvalue*",

where

  * "*tagname*" is a string

  * "*tagvalue*" is a single numerical or string item or a list of items between "(" and ")"
  parenthesis and separated by commas. e.g. (*cathode, anode*)

Format is free for the parameter assignements, provided that they are separated by
spaces. Everything which follows a '#' is considered as a comment and is disregarded.

For example::
  driftdiffusion
    {
     coupling = poisson
     nonlin_max_it = 70
     nonlin_rel_tol = 1e-10
     ls_max_step = 2
     ksp_type = bcgs
    }

Here and in the whole input file a string item can include a combination of characters,
special characters and numbers, but not spaces; if a space is found , the string item is
taken as terminated.

.. note:: |more|

   The input file is composed by the following sections:

      **Device, Scale, Models, Physics, Solver, Simulation**

   which will be described in the following.

Device section
--------------

  $Device
    {
     Region buffer
       {
        .........
       }
     Region barrier_1
       {
        ............
       }
     .................
    }
    
In "**Device**" section, two kinds of block can be present: the **Region** block and the 
**Cluster** block.

The **Region** blocks contain the description of the device in continuous media approach; 
the **Cluster** blocks define each a group of regions (mesh_regions) even with
different physical properties, but to be treated together somewhere in the simulation
(e.g. quantum calculation). In this way it is possible to refer to the set of these regions
simply by the **Cluster** name.

Each **Region** block must be preceded by the keyword "**Region**", followed by the
(single-word) name of the **TiberCAD Region**. The name of the **TiberCAD Region** 
can coincide with the name of a *mesh_region*, as defined during the modeling of the
device; in this case, if the keyword *mesh_regions* is absent, the **TiberCAD Region** will
be associated to the mesh region identified by the name assigned to the  **TiberCAD
Region**. Otherwise, the **TiberCAD Region** will be associated to the mesh regions
specified by the keyword *mesh_regions*::

  Region QWell
    {
     mesh_regions = (4,5)
     # mesh_regions = 4
     structure = wz
     y-growth-direction = (1,0,-1,0)
     z-growth-direction = (-1,2,-1,0)
     x-growth-direction = (0,0,0,1)
     material = GaN
     doping = 1e17 doping_type = donor doping_level = 0.025
    }

Here are the description of the available keywords for a **Region** block.

* *material* (mandatory): name of the material associated to the present region, e.g **Si**;
it may be a ternary alloy, e.g **AlGaAs**, in this case keyword *x* described in the following
has to be present.

* *x* : alloy concentration, expressed as the molar fraction of the first component of
the alloy; e.g. to express an alloy :math:`Al_xGa_{1-x}As` with molar fraction x = 0.2, that is
:math:`Al_{0.2}Ga_{0.8}As`, we select **AlGaAs** for the keyword material, and 0.2 for the keyword x,
thus we write x = **0.2**.

* *mesh_regions* : (a list of) region ID(s) or physical name(s) as specified in the meshing
program.

* *structure* : crystal structure (wz = wurtzite, zb = zincblend)

* *x-growth-direction, y-growth-direction, z-growth-direction* : Bravais vectors with Miller
indexes for wurtzite crystal (4 element vectors) or zincblende crystal (3 element vectors).

* *doping* : doping concentration [:math:`cm10^{-3}`]

* *doping_type* : donor or acceptor

* *doping_level* : energy level of the dopant [eV]

Each **Cluster** block must be preceded by the keyword "**Cluster**", followed by the
(single-word) name of the **Cluster**::

  Cluster Quantum_1
    {
     mesh_regions = (3,4,5)
    }
  
* *mesh_regions* (mandatory): list of physical regions (region IDs or physical names) or
TIBERCAD region names, as specified in the meshing program, to be grouped in the
cluster.

**Regions** and **Clusters** represent the macroscopical description of the device or structure 
to be be simulated in **TiberCAD**. In the rest of the input file, the physical regions
associated to **Models** or **Atomistic** descriptions will be indicated by means of the
**TiberCAD_Region** and **Cluster** names.

Scale section
-------------

The section "**Scale**" is dedicated to the optional definition of Non-Continuous Media
regions for the device: these regions will be described and studied at a different scale (e.g.
atomistic, circuit level lumped model, etc.) As for now, just the atomistic description
is implemented. **Atomistic** blocks, if present, specify a possible atomistic description
associated to one or to a group of physical regions described by **Region**  blocks. So,
for each **Atomistic** block defined in Input file, an atomic structure description will
be generated and used to solve a simulation problem with an atomistic approach. The
association to the physical (macroscopic) regions of the device allows the implementation
of multi-scale calculations.

Each **Atomistic** block must be preceded by the keyword "**Atomistic**", followed by
the (single-word) name of the atomistic region::

  Atomistic TB_1
    {
     physical_regions = (barrier_1, qwell, barrier_2)
     # physical regions to be described
     # with atomistic model.
    }
  
Here are the description of the available keywords for an **Atomistic** block.

* *physical_regions* (mandatory): list of the physical regions (**TiberCAD Regions** or 
**Clusters**) of the device associated to an atomistic description. all (default) is used to
specify all the physical regions.

* *path* (optional): path for importing an atomistic structure from an external file. xyz
and gen formats are supported, and are automatically recognized by file extension. Each
of the atom positions is imported as is, so the atom coordinates must be consistent with
the geometry of the device.
If no path is specified, the **TiberCAD Atomistic Generator** builds the atomistic
structure; it is constructed as a bulk crystal structure, covering with proper atomic
species the physical regions and taking in account the dimension of the problem (up to
now 1D structures are supported)............

Atomistic Generator options, to be put in the Atomistic section, are described in the
following.

* *reference_region* (mandatory): the Atomistic Generator can only build pseudomor
phical heterostructures. A reference region must be defined to specify from which region
(**TiberCAD Regions**) to get structure parameters such as lattice constants, which
depend on the material defined in the reference region.

* *passivation* (optional): no is default option. If set to *no*, no passivation is performed.
If *yes* is specified, a hydrogenation of the structure is performed, taking into account the
structure periodicity. Up to now, hydrogenation is supported fo Silicon structures.

* *preserve* (optional): Default is *none*. In some cases it is requested to build a structure
in which the atom basis or the conventional cell has be preserved, regardless to mesh
geometry. If *none* is specified, no conservation is performed and only atoms strictly belonging 
to geometrical regions are put in the atomic structure. If *lattice* is specified, atom 
basis is preserved (e.g. to preserve anion/cation couples). If *conventional* is specified,
conventional cell is preserved. ..............

* *y_lenght* (optional): Atomistic Generator builds the minimum periodical structure
along *y* and *z* directions. If y lenght is specified, the structure will be at least y lenght
sized along *y_growth_direction*. Exact lenght is internally defined in order to keep struc-
ture periodicity.

* *z_lenght* (optional): same as above, for the z direction.

Models section
--------------

::

  $Models
    {
     model driftdiffusion
      {
       ..........
       BC_Regions
         {
          BC_Region cathode
            {
             ............
            }
          BC_Region anode
            {
             ..............
            }
         }
     }
    model macrostrain
    {
     ..............
  
In "**Models**" section, one or more model-blocks must be present: each model-block
must be preceded by the keyword "**model**", followed by the (single-word) model name.
This must be the name of one of the TIBERCAD simulation models.

Here are the simulation models implemented until now:

* **driftdiffusion**: Poisson-driftdiffusion transport of electrons and holes

* **thermal**: Heat balance simulation

* **excitontransport**: Exciton transport model

* **macrostrain**: Calculation of Elastic deformations in heterostructures

* **efaschroedinger**: Envelop Function Approximation (EFA) solution of single 
 particle Schroedinger equation for electrons and holes

* **quantumdensity**: Calculation of quantum density of electrons and holes.

* **quantumdispersion**: Dispersion of quantized states in k space

* **opticskp**: Optical properties (optical kp matrix elements)

* **opticalspectrum**: Emission spectrum (with k-space integration)

For a complete description of these simulation models, see the next chapters.

Each model-block can contain some optional blocks, to be written in the following
order:

* one "options" block, preceded by the keyword "**options**". This block can contain 
general options for the present model.

* one or more **physical_model** blocks: each physical model block must be preceded
by the keyword "**physical_model**", followed by the (single-word) name of the physical model. 
Each physical_model block can contain parameters relevant to a specific
model of a physical property or quantity related to the present model.

* one or more Boundary Condition regions blocks (**BC_regions-block**). The BC_regions-block
 must be preceded by the keyword "**BC_Regions**" and it is composed by one or more parameters-blocks,
 each preceded by the keyword "BC_Region" followed by the (single-word) name of the 
 boundary condition region. This parameters-block can contain the possible description of the boundary region.

These optional blocks must be strictly in this order: first the **options**, then the
**physical_model**, and finally the **BC_regions-blocks**. 

A detailed description of the possible parameters for these blocks follows.

Options block
^^^^^^^^^^^^^

::

  options
    {
     simulation_name = dd1
     # physical_regions = all
     physical_regions = (channel_1, channel_2)
    }
  
* *simulation_name*: user-defined name of the particular instance of the simulation
model defined for this block. More than one simulation (with different name and properties) 
can be defined, in separated model blocks, which refer to a same TIBERCAD
simulation model. If simulation name is not assigned, by default the TIBERCAD modelname 
is taken as current simulation name.

* *physical regions*: (list of) physical region(s) to which the present simulation model
will be applied. Physical region(s) are the **TiberCAD regions** or **clusters** as defined
in **Device** section
Default value is "all" (all physical regions of the device). In a list, the names must
be separated by comma and enclosed between '(' and ')' parenthesis

Physical model block
^^^^^^^^^^^^^^^^^^^^

::

  physical_model recombination
    {
     model = SRH
    }
  
The following options can be applied to any physical_model:

* **name** a user defined name for the model for unique identification

* **restrict_to_region** a single region name or a list of regions where the model should be used. 
This allows to use different implementation of the same type of model in
different regions.

.. _BC_region:

BC_region block
^^^^^^^^^^^^^^^

::

  BC_Region anode
    {
     BC_reg_numb = 2
     type = ohmic
     # voltage = 0.0
     voltage = @Vb[1.5]
    }
  
* *BC_region*: name of the present boundary region; it can be the name of a boundary
(physical) region specified in the meshing program (GMSH or ISE-TCAD).
* *BC_reg_numb*: BC_region ID(s) as specified in the meshing program (GMSH). If this 
keyword is not present, it is assumed that the mesh region associated to this **TiberCAD BC_region** is given by the name assigned to *BC_region*.

* *type*: type of boundary condition: ohmic, schottky,... , substrate (for strain calculations).

* *voltage*: value of voltage [V] applied to the present BC region (for ohmic and schottky
BCs); it can be the value of a sweep variable 'Vb', indicated with @Vb. A possible default
value can be indicated in parenthesis: @Vb[1.5]
* *zero_grad_fermi_h, zero_grad_fermi_e*: if true set Neumann b.c. to the fermi level in
the b.c. region.

If type is substrate (for **strain** calculations):

* *material* : name of material in the substrate region.
* *structure* : crystal structure ( wz = wurtzite, zb = zincblend)
* *x-growth-direction, y-growth-direction, z-growth-direction*: Bravais vectors with Miller 
indexes for wurtzite.

Solver section
--------------

----
::

  $Solver
    {
     driftdiffusion
       {
        nonlinear_solver = tiber
        ksp_type = bcgsl
        nonlin_rel_tol = 1e-12
        nonlin_abs_tol = 1e-15
        nonlin_step_tol = 1e-2
        # nonlin_rel_tol = 1e-12
        lin_rel_tol = 1e-6
        nonlin_max_it = 30
        ls_max_step = 1
        # pc_type = lu
        pc_type = composite
       }
    ...................

----
 
In this section one may choose and define the setting parameters for the numerical 
solvers to be applied to the specified simulations: the section is organized in blocks, each 
one preceded by a block-name. This block-name may be

#. the name of one of the user-defined simulations: in this case, the solver parameters 
defined in each block will refer only to that particular simulation.

#. the name of a TIBERCAD **Model**: only in this case, the settings will be applied to 
*all* the simulations belonging to that **Model**

For details on the available parameters for each **Model**, see the relevant chapter in this
Guide.

In addition, two special (optional) blocks may be present: the **Sweep** block and the
 **Selfconsistent** block

* **Sweep** block, preceded by the keyword **"Sweep"**. This block may contain one 
or more subblocks, each one defining a set of calculations applied to a boundary
region (e.g. a set of bias values to be assigned to a drain contact of a MOSFET
for the calculation of an output drain IV characteristic), in this Guide referred to
as *sweep* calculation.
Each *sweep* definition must be preceded by its user-defined name (e.g. *sweep_1*).
( see :ref:`Listing1` )
In the case of a single *sweep* subblock, for backward compatibility, it is also allowed
to define the *sweep* feature by means of the keyword *sweep* (lower case), in the
following way::

  sweep
    {
     simulation = driftdiffusion
     variable = Vb
     start = 0.0
     stop = 4.0
     steps = 80
     plot_data = true
     plotvariable = current
    }

The following keywords are defined for this feature:

**variable**: name of the variable to which the *sweep* is applied: its value is assigned
to a quantity (e.g. voltage) in a BC_Region section to perform the *sweep* calculation
( see :ref:`BC_region` ).
**start, stop, steps**: sweep starts from *start* value, is repeated *steps*  times and stops 
in *stop*
**simulation**: name of the simulation (model) associated to the sweep calculation;
it may be the name of another sweep defined in the same block.
**plotvariable** (obsolete): specify the integrated quantity to be calculated during the
sweep and that will be shown in the output file  *sweep_modelname_sweepvariable.dat* ,
eg. ``sweep_driftdiffusion_Vb.dat`` for a sweep of current calculation on the variable
Vb (typically a contact voltage).

.. _Listing1: Listing1

----

::

  Sweep
    {
     sweep_1
      {
       simulation = driftdiffusion
       variable = Vd
       start = 0.0
       stop = 2.0 #0.1
       steps = 20 #50
       # plot_data = true
      }
  
     sweep_2
      {
       variable = Vg
       start = -0.1
       stop = 1.0 #0.5
       steps = 11 # 6
       simulation = sweep_1
       # simulation = driftdiffusion
       plot_data = true
      }
    }
    
----

Listing 1: Example of Sweep section

**plot_data**: default is *false*; if it is set to true, then output data will be written for
each step of the sweep calculation, otherwise just the results for the final step will 
be present in the output.
Once a *sweep* calculation has been defined, it is treted as a special case of simulation 
and may be executed as an usual *simulation*: by adding it in the **solve** list, e.g.

  *solve* = *sweep_drain* ( see :ref:`SimulationSection` ).

* **Selfconsistent block**, preceded by the keyword **"Selfconsistent"** . This block may
contain one or more subblocks, each one defining a self-consistent calculation based
on two different simulation models (e.g. driftdiffusion and excitontransport). The
definition of a *selfconsistent* calculation must be preceded by its user defined name
(e.g. *converse_piezo*):

.. _Listing2: Listing 2

----

::

  Selfconsistent
    {
     converse_piezo
       {
        flavour = relaxation
        simulations = (driftdiffusion, strain)
       }
     dd_excitons
       {
        flavour = relaxation
        simulations = (driftdiffusion, exciton)
       }
    }

----
 
Listing 2: Example of Selfconsistent section

The following keywords are defined for this feature:

**simulations**: the list of simulations to be performed self-consistently.

**flavour**: specifies **broyden** or **relaxation** approach

Once a *selfconsistent* calculation has been defined, it may be executed as an usual
simulation: by adding it in the **solve** list, e.g. *solve = dd_excitons* ( see :ref:`SimulationSection` ), or
even in the **sweep** section, with, for example, *simulation =converse_piezo* . In 
both cases, the specified simulations will be executed in a self consistent way.

In the case of a single self-consistent calculation block, for backward compatibility,
it is also possible to define the *self-consistent* feature by means of the keyword 
**selfconsistent** (lower case), in the following way::

  selfconsistent
    {
     flavour = relaxation
     simulations = (driftdiffusion, excitons)
    }


Physics section
---------------

In this section several physical parameters can be entered, in addition to or overwriting 
the material parameters present in the material files. The section is organized in blocks, 
each one preceded by a block-name.
   
   This block-name may be:

#. the name of one of the user-defined *simulations*: in this case, the physical param-
eters defined in each block will refer only to that particular simulation.

#. the name of a TIBERCAD **Model**: only in this case, the settings will be applied to
ALL the *simulations* belonging to that **Model**

For details on the available parameters for each **Model** , see the relevant chapter in this Guide.

  For example in this case:
  **strain_simulation** is used to specify the simulation that provides strain data, in case 
of strained systems ( see :ref:`PhysicSection` ).

.. _SimulationSection:

Simulation section
-------------------

In this section one can specify several general parameters and settings for the actual 
calculation to be run, such as the mesh file to be used, the dimension of simulation, the 
process-flow of simulation, etc.

* *searchpath*: path for material files:
----

::

  driftdiffusion
    {
     strain_simulation = strain
    }

----

* *meshfile*: name of mesh file. |warn| N.B.: the extension is mandatory! ( *.grd* for ISE-TCAD, 

*.msh* for GMSH mesh file v.1 and v.2.0 )

* *mesh_units*: units of measurements used in the meshing (relative to meters): e.g.,
:math:`10^{-6}` for :math:`\mu m`

* dimension: dimension of simulation (1,2,3)

* temperature: temperature of the system [K]

* solve: list of simulations to be executed, in the order of execution; if the list contains 
**"sweep"** , a sweep is performed as specified in *sweep* block in the *Solver* section::

  solve = (strain,driftdiffusion, quantum_electrons, quantum_holes)

* resultpath: path for output directory

* output_format : format of the output data: gmv for GMV, ise for **Tecplot**, *grace* for
**xmgr** (ascii data column type), *vtk* for **Paraview**.

* plot : list of output variables which are calculated and available in output files. See 
the corresponding chapters for the available output variables for each model.

Output description
------------------

At the end of the execution, the program will write the results of the simulation in the
directory specified by *resultpath*, with the format specified by *output_format*.
 The output variables are specified in the list plot.

TiberCAD output is divided in three classes: **nodal, elemental** and **integrated** quantities.

**Nodal** quantities are all the quantities associated with the nodes of the mesh, such
as Fermi level, electron and hole density, conduction and valence band, etc. The output
values for these quantities are reported in the files *modelname_nodal.ext*, where *modelname* 
is the simulation model used for the calculations and ext is the extension of the chosen file format.

In the case a sweep calculation is performed and the **plot_data** keyword is set to *true*,
the output files are of the kind *modelname_nodal_sweepvariable_step.ext*, where *sweep variable* 
is the variable with respect to which the sweep is performed (e.g. gate voltage)
and *step* is the value of this variable at that step; e.g *driftdiffusion_nodal_Vb_0.000.dat*
for the result at the step Vb = 0.0.

**Elemental** quantities are all the quantities associated with the elements of the mesh,
such as current density. The output values for these quantities are reported in the files
:file:`modelname _elemental.ext`.
In the case a sweep calculation is performed and the **plot_data** keyword is set to *true*, 
the output files are of the :file:`kind_modelname_elemental_sweepvariable_step.ext`. e.g.
:file:`driftdiffusion_elemental_Vb_1.150.dat`.

**Integrated** quantities are the quantities which are not associated to the mesh but are
obtained by an integration on real or reciprocal space, for example current at the contacts
of a diode or quantized energy levels in a quantum well. These **Integrated** quantities
are displayed in separated files, with the format simname.ext, e.g :file:`quantum_electrons.dat`,
where simname is the name of the model (simulation) associated to the results. If a
sweep is performed, the output file gets the :file:`format_sweep_simname_varname.ext`, where
varname is the variable with respect to which the sweep is performed, for example
:file:`sweep_driftdiffusion_Vb.dat`. 

Inside the file, output values for all the steps of calculation are shown.

Finally, a last class of output files is the **Materials** output. These files contain
the information about the physical regions of the device; for each class of simulation, a
different material file is produced, containing all and only the mesh regions associated
to that simulation model. The file has the :file:`format_simulationname_materials.ext`, e.g.

:file:`driftdiffusion_materials.dat`


Example of Input file
---------------------

Here is an example of the input file template::

  # Description of the device physical regions
  $Device
    {
     # Syntax:
     #
     # Region "Tiber_region"42 CHAPTER 4. INPUT FOR TIBERCAD
     # mesh_regions = "list gmsh region ID/names" | "list ISE_TCAD region names"
     #
     # if mesh_regions is empty -> mesh_regions = "Tiber_region"
     #
     Region buffer
       {
        mesh_regions = 1
        structure = wz
        y-growth-direction = (1,0,-1,0)
        z-growth-direction = (-1,2,-1,0)
        x-growth-direction = (0,0,0,1)
        material = GaN
        doping = 1e15
        doping_type = donor
        # doping_level = 0.025
       }
     Region barrier_1
       {
        mesh_regions = (2,3)
        structure = wz
        y-growth-direction = (1,0,-1,0)
        z-growth-direction = (-1,2,-1,0)
        x-growth-direction = (0,0,0,1)
        material = AlInN
        x = 0.80 #
        doping = 1e15
        doping_type = donor
        # doping_level = 0.025
       }
     Region QWell
       {
        mesh_regions = 4
        structure = wz
        y-growth-direction = (1,0,-1,0)
        z-growth-direction = (-1,2,-1,0)
        x-growth-direction = (0,0,0,1)
        material = GaN
        doping = 1e15 #
        doping_type = donor #
        # doping_level = 0.025
       }
     Region barrier_2
       {
        mesh_regions = (5,6)
        structure = wz
        y-growth-direction = (1,0,-1,0)
        z-growth-direction = (-1,2,-1,0)
        x-growth-direction = (0,0,0,1)
        material = AlInN
        x = 0.80 #
        doping = 1e15 #
        doping_type = donor
        # doping_level = 0.025
       }
     # Cluster = group of mesh_regions with DIFFERENT material (in general)
     # Syntax:
     #
     # Cluster "Tiber_cluster"
     # mesh_regions = "list gmsh region ID/names" | "list ISE_TCAD region names"
     #
     # "Tiber_cluster" to be used in Models section
     Cluster Quantum_1
      {
       mesh_regions = (3,4, 5)
      }
    }
  #Definition of the description scale (only for not-Continuous Media regions )
  # Syntax:
  #
  # "level" "scale_cluster"
  # "level" = "Atomistic | ..."
  # physical_regions = "list (Tiber_region | Tiber_cluster) "
  #
  #******* $Scale section is optional *********
  #
  $Scale
    {
     Atomistic TB_1
      {
       physical_regions = (barrier_1 , QWell , barrier_2 )
      }
     Atomistic TB_2
      {
       physical_regions = ...........
      }
    }
  # Definition of Simulation Models and associated Boundary Conditions
  $Models
    {
     model driftdiffusion
      {
       options
        {
         simulation_name = dd1
         physical_regions = all
        }
       physical_model recombination
        {
         model = srh
        }
     physical_model recombination
      {
       model = direct
       C = 1.1e-8
      }
     BC_Regions
       {
        BC_Region cathode
          {
           BC_reg_numb = 1
           type = ohmic
           voltage = 0.0
          }
        BC_Region anode
          {
           BC_reg_numb = 2
           type = ohmic
           voltage = 0.0
          }
       }
     }
    model macrostrain
      {
       options
        {
         simulation_name = strain
         physical_regions = all
        }
     BC_Regions
       {
        BC_Region substr
         {
          BC_reg_numb = 1
          type = substrate
          material = GaN
          structure = wz
          y-growth-direction = (1,0,-1,0)
          z-growth-direction = (-1,2,-1,0)
          x-growth-direction = (0,0,0,1)
         }
       }
     }
   model efaschroedinger
     {
      options
        {
         simulation_name = quantum_electrons
         # Syntax:
         #
         # physical_regions = "list (Tiber_region | Tiber_cluster) "
         physical_regions = Quantum_1
        }
     }
   model efaschroedinger
     {
      options
        {
         simulation_name = quantum_holes
         physical_regions = Quantum_1
        }
      }
   }
  # Definition of Model-dependent Solver parameters
  $Solver
    {
     driftdiffusion
       {
        coupling = poisson
        ksp_type = bcgsl
        nonlin_abs_tol = 1e-10
        nonlin_step_tol = 1e-2
        # nonlin_rel_tol = 1e-12
        lin_rel_tol = 1e-6
        nonlin_max_it = 30
        # local_density_scaling = true
        # ls_type = none
        discretization = fem
        ls_max_step = 1
        # pc_type = lu
        pc_type = composite
        integration_order = 2
        # relaxation_factor = 0.5
       }
     macrostrain
       {
        substrate = substr
       }
     efaschroedinger
       {
        x-periodicity = false
        Dirichlet_bc_everywhere = true
        # particle = hl
        number_of_eigenstates = 30
        # model = conduction_band #eff mass cb
        poisson_model_name = driftdiffusion # potential from driftdiffusion
        strain_model_name = macrostrain
        convergent_density = true4.9. EXAMPLE OF INPUT FILE 49
       }
     quantum_electrons
       {
        particle = el
       }
     quantum_holes
       {
        particle = hl
       }
     }
  # Definition of Model dependent physical parameters
  $Physics
    {
     driftdiffusion
      {
       statistics = FD
       strain_simulation = macrostrain # default driftdiffusion model including
       # local strain obtained from "macrostrain"
      }
     quantum_electrons
       {
        particle = el
        model = conduction_band #eff mass cb
       }
     quantum_holes
       {
        particle = hl
        model = kp # k.p for valence band
        kp_model = 6x6
       }
    }
  # Definition of model-indipendent parameters of the Simulation
  $Simulation 
    {
     searchpath = ../../materials
     mesh_units = 1e-9 #nm !!
     dimension = 1
     meshfile = test.msh
     temperature = 300
     solve = (strain, driftdiffusion, quantum_electrons,quantum_holes )
     resultpath = output_
     plot = (Ec, Ev, QFermi_e, QFermi_h, EField ,
     eDensity, hDensity, eCurrent, hCurrent,
     Current, NetRecombination, eMob, hMob, T,
     strain, polarization, xEffPot, xDensity,
     xMob, ExcitonRecombination, EigenFunctions,
     EigenEnergy,EnergyLevels, xCurrent)
     output_format = grace
    }



.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes