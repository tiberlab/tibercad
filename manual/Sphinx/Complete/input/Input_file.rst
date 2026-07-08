


.. _InputFileGetting:




tiberCAD Input File  
=================================================


Input for tiberCAD is composed by an input file e.g. ``input.tib`` and a mesh file generated 
by a mesher software: as for now, mesh files from ``GMSH`` (\*.msh)
and from Synopsys TCAD (\*.grd) are supported.
Be sure that the material files are in the correct directory.


To run the program, type: 

  ``tibercad`` *input_file_name*

Description of Input file structure
-------------------------------------------

A valid input file for tiberCAD is a text file with the structure described in the following.
In the whole input file, everything following a "#" is considered as a comment and is
ignored; blank lines can be present anywhere and are ignored, too.
An input file is composed by several blocks :

A block is enclosed between curly braces ("{", "}") and may include one or more blocks.
Each block has a header made of one or two keywords.
Each block may contain zero or any number of parameter assignments of the
form ``tagname = tagvalue``, where

  * ``tagname`` is a string

  * ``tagvalue`` is a single numerical or string item or a list of items in parentheses
    and separated by commas. e.g. ``(cathode, anode)``

Format is free for the parameter assignments, provided that they are separated by
spaces. Everything which follows a "#" is considered as a comment and is ignored.

For example::

  electron_mobility field_dependent
  {
    region = buffer
    # type = doping_dependent
    low_field_model = constant
  }

Here and in the whole input file a string item can include a combination of characters,
special characters and numbers, except spaces; if a space is found, the string item is
assumed to be terminated.

The input file is composed by the following main classes of blocks:

    ``Device``, ``Module``, ``Simulation``

which will be described in the following.

Device section
-------------------

::

  Device hemt
  {
    meshfile = hemt.msh

    Region buffer
    {
      ...
      Doping
      {
        ...
      }
    }

    Region barrier
    {
      ...
    }
    ...
  }

The ``Device`` section includes the geometrical description of the device to be simulated.
An optional device name can be associated to the device object after the ``Device`` keyword.
In the ``Device`` section, two types of blocks used for the description of the device geometry can be present: ``Region`` and 
``Cluster`` blocks. 

Outside of these blocks, general options common to all the device can
be defined. The most important one is the specification of the mesh file, which is mandatory.



 ``meshfile`` : string
     name of the mesh file, including file name extension
     (``*.grd`` for Synopsys devise, ``*.msh`` for ``GMSH`` mesh files)

 ``mesh_units`` : double
     units used in the meshfile in meters, default is ``1e-6`` corresponding to micrometers


.. warning::
             For a *1D* simulation the geometrical model (mesh) has to be drawn on the *x axis*.
             For a *2D* simulation the geometrical model (mesh) has to be drawn in the *xy-plane (z=0)*. 
             Any other orientation will produce wrong results.



The ``Region`` blocks contain the description of the device in continuous media approach.
The ``Cluster`` blocks define logical groups of regions, which may have different materials or
different physical properties. In this way it is possible to easily refer to sets of regions
by using the cluster name.

Two further blocks may be present in the ``Device`` section:
``Atomistic`` and 
``Parallel`` blocks. 

The ``Atomistic`` block contain the description  of an atomistic structure which can be related to one or more Regions of the Device.
The ``Parallel`` block defines the parameters which control a parallel execution of the Device simulation. 


Region block
^^^^^^^^^^^^

``Region`` blocks are started with the keyword ``Region`` , followed by the
name of the tiberCAD region. The name of the tiberCAD region 
can coincide with the name of a mesh region, as defined during the modeling of the
device. In this case, if the keyword ``mesh_regions`` is not used, the tiberCAD region will
be associated to the mesh region identified by the given name.
Otherwise, the tiberCAD region will be associated with the mesh regions
specified using the keyword ``mesh_regions``.
For example, in  the  following, the tiberCAD ``Region``  *QuantumWell* will comprise the  two  mesh regions *well1* and *well2*, defined during  the  geometrical modeling  of  the  device  

::

  Region QuantumWell
  {
    mesh_regions = (well1, well2)

    Doping
    {
      density  = 1e17
      type = donor
      level  = 0.025
    }
  }

In  the following, instead, the  ``Region``  *well1* takes  its  name  from the mesh  region *well1*  it  refers to

::

  Region well1
  {
    
    Doping
    {
      density  = 1e17
      type = donor
      level  = 0.025
    }
  }

 
The available keywords inside a ``Region`` block are the following:

  ``material`` : string
             name of the material associated to the region, e.g ``Si`` or ``AlGaAs``

  ``x`` : double, 0 < ``x`` < 1
          alloy concentration, expressed as the molar fraction of the first component of the alloy.
          To create e.g. an alloy :math:`\mathrm{Al}_{0.2}\mathrm{Ga}_{0.8}\mathrm{As}`,
          we specify ``AlGaAs`` for the keyword ``material``, and ``0.2`` for ``x``

  ``mesh_regions`` : string, list
          a list of region names as specified in the meshing program.

  ``[x,y,z]-growth-direction``  : 3-/4-tuple
          Bravais vectors with Miller or Miller-Bravais indices indicating the crystal directions
          along the x, y and z axes. One may specify only two or even a single direction, as the
          others will be computed internally. Based on the type of braces, crystal directions or
          lattice planes will be used, based on standard notation using ``[]`` for directions and
          ``()`` for planes.

  ``euler_angles`` : 3-vector
          As an alternative to Miller indices, the crystal orientation can be specified via Euler
          angles describing the rotation from the standard crystal orientation to the one in the
          calculation coordinate frame. We use the convention :math:`Z_\gamma Y_\beta Z_\alpha`,
          and thhe order of the input arguments is :math:`(\alpha, \beta, \gamma)`.

.. ``structure`` : string
..          crystal structure (wz = wurtzite, zb = zincblend)


.. tip:: A common material, crystal structure and growth directions can be defined for all
          device regions by defining them outside of the ``Region`` blocks.



The optional subblock ``Doping`` as in the example above contains the keywords:

  ``type`` : string
         The dopant type. Can be ``donor`` or ``acceptor``.

  ``density`` : double 
         The doping concentration in cm\ :sup:`-3`.

  ``level`` : double
         The energy level of the dopant given as the distance from the conduction band edge (for donors)
         or from the valence band edge (for acceptors) in eV.

  ``alpha`` : double
         The parameter for the doping density dependence of the ionization energy, given as
         :math:`\Delta E_d = \Delta E_{d,0} - \alpha N_d^{1/3}`, where :math:`\Delta E_d` corresponds to the value
         of ``level``, and :math`N_d` is the density given by ``density``. The default value is 0.0.


  ``g`` : integer
         Level multiplicity. Defaults to 2 for donors and 4 for acceptors.


Cluster block
^^^^^^^^^^^^^

The definition of cluster blocks must be preceded by the keyword ``Cluster`` , followed by the
name of the Cluster. For example

::

  Cluster Quantum
  {
    regions = (well, barrier1, barrier2)
  }

groups the regions ``well``, ``barrier1`` and ``barrier2`` in a logical unit with name ``Quantum``.
The  ``regions`` option is mandatory and contains the list of mesh or tiberCAD regions to be grouped
together.

Regions and clusters represent the macroscopical description of the device or structure 
to be be simulated in tiberCAD. In the rest of the input file, the regions
associated to the Modules will be indicated by means of tiberCAD region names and cluster names.

Atomistic structure
^^^^^^^^^^^^^^^^^^^

The definition of an atomistic structure can be done using the keyword ``Atomistic``. 

::
   
  Atomistic qw
  {
    regions = (all)
    reference_region = barrier1
  }

The atomistic generator produces an atomistic model over all the specified regions. Currently the structure
can only be a pseudomorphic structure. The reference region is used to define the crystal lattice for the 
whole structure. More details about the atomistic generator can be found in  :ref:`Atomisticgen`.


Parallel block
^^^^^^^^^^^^

The definition of the options for a parallel execution of the simulation can be done using the keyword ``Parallel``. 

::
   
  Parallel
  {
    mpi_processes_per_device = 2 # 4 
    # let's do FEM calculations in serial
    mpi_processes_per_mesh   = 1
  }

More details  can be found in :ref:`Parallel`.





Modules
------------

::

  Module driftdiffusion
  {
    name = dd
    #regions = all

    plot = (Ec, Ev, Eg, eQFermi, hQFermi,
            eDensity, hDensity, eMobility, hMobility,
            Polarization, ElPotential, ElField,
            eCurrentDensity, hCurrentDensity,
            CurrentDensity, ContactCurrents,
            eConductivity, hConductivity)

    Solver
    {
      max_iterations = 15
      relative_tolerance = 1e-12
      step_tolerance
    }

    Physics
    {
      
  
      mobility
      {
        type = field_dependent
        low_field_model = doping_dependent
      }
     recombination (srh, auger) {}

    }

    Contact cathode
    {
      voltage = $Vd
    }

    Contact anode
    {
      voltage = 0.0
    }
  }

One or more module-blocks may be present: each module-block must be preceded by
the keyword ``Module``, followed by the (single-word) module name. This must be the
name of one of the tiberCAD modules.
Here are the Modules implemented until now:

 ``driftdiffusion`` : 
       Poisson-driftdiffusion transport of electrons and holes

 ``thermal`` : 
       Heat balance simulation

 ``elasticity`` : 
       Calculation of elastic deformations in heterostructures

 ``efaschroedinger`` : 
       Envelop Function Approximation (EFA) solution of single particle Schroedinger equation for electrons and holes 

 ``opticskp`` : 
       Optical properties (spontaneous emission  spectrum)

 ``DSC`` : 
       Simulation of a DSC solar cell

 ``vff`` : 
       Relaxation of an atomistic structure with valence force field 

 ``empirical_tb`` : 
       Atomistic quantum calculations with an Empirical Tight-binding solver

 ``sweep`` : 
       Parameterized execution of a module simulation (e.g. for the calculation
       of output current characteristics)

 ``selfconsistent`` : 
       coupled calculations of different simulation modules


Each module-block usually contains a list of general options, such as  ``plot``  and others
specific to each module. Then, two main blocks define the ``Physics`` and the ``Solver``
models and parameters for this module.

  **Solver** contains the solver parameters; depending on the Module, it can contain a
  LinearSolver definition subblock.

  **Physics** usually contains the definition of the *physical models* used in the simulation.


Examples of *physical models* are, for the **driftdiffusion**  module, **recombination**, **electron mobility**, **trap**, **polarization** 
and so on. A particular model is the **Boundary** model, which has an alias **Contact**
for driftdiffusion module, as  we will  see in  the  following.
The declaration of these models obey to the following syntax:

    *model_keyword* type_specifier
    
        <block>

where *model_keyword* is the name of the physical model to be declared (e.g. the *trap*
model), **type_specifier** is the name of a particular one among the available descriptions 
for that model and the  following block contains the options for the model. For example, for a *trap* model of **type** *acceptor* ::

  trap acceptor
    {
     region = buffer
     Nt = 7e16
     Et = 0.5
     reference = cb
    }

An alternative multiple declaration is possible if no options, other than default, are
declared:

    *model_keyword* (type specifier1, type specifier2,...)

  ``recombination (srh, direct, auger) { }``

In this example, several recombination models are defined (srh, auger, direct) each one
with default parameters.
For a detailed description of the models, please refer to the reference guide.

Two special modules are  the ``sweep`` and the ``selfconsistent`` modules

Module sweep
^^^^^^^^^^^^^^^^^^^^^^^^

::

  Module sweep
  {
    name = sweep_drain
    solve = driftdiffusion
    variable = $Vd
    start = 0.0
    stop = 2.0
    steps = 20
    plot_data = true
  }



Each ``sweep`` Module defines a set of calculations applied to a **boundary** region (e.g. a
set of bias values to be assigned to a drain contact of a MOSFET for the calculation of
an output drain IV characteristic), in this Guide referred to as ``sweep`` calculation.

The following keywords are defined for this feature:

  ``variable`` 
     name of the variable to which the ``sweep`` is applied: 

Writing e.g. ``variable = $Vg`` means that a ``sweep`` on the variable ``$Vg`` will be performed.
(see  :ref:`Variable`).

Further options controlling a ``sweep`` are:

  ``name``
     A user defined name. This is important if several ``sweeps`` are defined.

  ``solve``  
     names of the simulations (modules) to be solved at each ``sweep`` step.
     The names of other ``sweeps`` may be provided to realize nested ``sweeps``.

  ``start, stop, steps`` 
     The ``sweep`` goes from ``start`` to ``stop`` in ``steps`` steps.

  ``min_step``
     The minimum allowed step size.

  ``max_step``
     The maximum allowed step size.

  ``plot_data``  
     If set to ``true`` or ``each``, output data will be written after
     each step of the ``sweep``. Other options are ``last`` or ``none``.
     Alternatively, one can specifiy explicitly the values for which to write results.

Once a ``sweep`` calculation has been defined, it is treated as an ordinary simulation
and may therefore be used just like any other simulation by adding it to a ``solve`` statement (see below).

Module selfconsistent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this Module it is possible to define a self-consistent calculation based on two different
simulation modules (e.g. ``driftdiffusion`` and ``thermal``).

::

  Module selfconsistent
  {
    name = sc_all
    solve = (tb, dens_el, dens_hl, dd)
    max_iterations = 5
    absolute_tolerance = 1e-3
    relative_tolerance = 1e-6
  }

In ``solve`` the list of simulations to be performed self-consistently is specified.

See  also :ref:`Selfcons`.   



Simulation section
--------------------------

::

  Simulation
  {
      temperature = 300
      verbose = 3
      solve = (strain, dd, sweep_g)
      logfile = output_Nt5e16_Al0.215_SiN/hemt.log
     
      resultpath = output_Nt5e16_Al0.215_SiN
      output_format = vtk
  }


In this block one can specify several general parameters  such as the temperature and the settings for the actual
calculation to be run,that is  the process-flow of the simulation.

  ``searchpath`` : 
     path for material files, default is  the system-defined material directory.

  ``temperature`` : 
     temperature of the system [K]

  ``solve`` : 
     list of simulations to be executed, in the order of execution

     ``solve = (strain,driftdiffusion, quantum_electrons, quantum_holes)``

  ``resultpath`` : 
     path for output directory

  ``verbose`` :
    a number to control the verbosity of output to the screen and the log file


Output description
-------------------------

At the end of the execution, the program will write the results of the simulation in the
directory specified by  ``resultpath``, with the format specified by ``output_format``. 
The output variables for each Module are specified in the list *plot*.

tiberCAD output is divided in two classes: **mesh-based** and **mesh-independent**
quantities.

The  available  formats are 

  * ``vtk`` for 2D and 3D  output data

  * ``grace`` for 1D output  data

Output results from 2D or 3D simulations may  be  visualized with  the open source visualization and post-processing  tool **paraview**: http://www.paraview.org

Output results from 1D simulations are  ascii  data  files which  can  be visualized e.g.  with the  open source plotting tools: 

**xmgrace** : http://plasma-gate.weizmann.ac.il/Grace
(for Linux)

**AptPlot** :
http://www.aptplot.com/aptplot/
(in java,  for  Win and  Linux)



Mesh-based quantities
^^^^^^^^^^^^

The ``mesh-based`` quantities are all the quantities associated with the nodes of the mesh,
such as Fermi level, electron and hole density, conduction and valence band, etc., together
with all the quantities associated with the elements of the mesh, such as current density.

The output values for these quantities are reported in the files *simname_msh.ext*,
where *simname* is the  Module used for the calculations and *ext* is the extension
of the chosen file format, e.g. **vtu** for paraview output.

  ``strain_msh.vtu``

In the case a ``sweep`` calculation is performed and the **plot_data** keyword is set to
*true*, the output files are of the kind *simname_sweepvariable_stepvalue_msh.ext*, where *sweepvariable* is the variable with respect to which the sweep is performed (e.g. gate voltage) and *stepvalue* is the value of this variable at that step; e.g the result at the
step ``Vb = 1.1`` will be found in the file:

  ``dd_Vb_1.1_msh.vtu``

Mesh-independent quantities
^^^^^^^^^^^^^^^^

The ``mesh-independent`` quantities are the quantities which are not associated to the
mesh, for example current at the contacts of a diode or quantized energy levels in a
quantum well. 

These mesh-independent quantities are displayed in separated files,
with the format ``simname.ext``, e.g ``quantum_electrons.dat`` , where ``simname`` is the name of
the Module associated to the results and ``ext`` is the filename extension (usually ``dat``). If a sweep is performed, the output file
gets the format ``sweep_name_simname.ext``,   
where ``sweep_name`` is the name of the ``sweep``
performed, for example

  ``sweep_drain_driftdiffusion.dat``
The mesh independent results for every sweep step are stored in this file.




Material database
----------------------------------







..  math::
    :nowrap:
    :label:

     \begin{table}[!ht]
     %\center
     \begin{minipage}[t]{8cm}
     \begin{tabular}{l|c|l}
     \multicolumn{3}{c}{\textbf{Material list}} \\
     \hline
     \textbf{Name}  & \textbf{Crystal structure} & \textbf{Type}  \\
     \hline
     \hline
     \texttt{Air} & - & - \\
     \texttt{AlAs} & zb & binary \\
     \texttt{AlAsP} & zb & ternary \\ 
     \texttt{AlAsSb} & zb & ternary  \\
     \texttt{AlGaAs} & zb & ternary  \\
     \texttt{AlGaN} & wz & ternary  \\
     \texttt{AlGaP} & zb & ternary  \\
     \texttt{AlGaSb} & zb & ternary  \\
     \texttt{AlInAs} & zb & ternary  \\
     \texttt{AlInGaAs} & zb & quaternary  \\
     \texttt{AlInGaP} & zb & quaternary  \\
     \texttt{AlInN} & wz & ternary  \\
     \texttt{AlInP} & zb & ternary  \\
     \texttt{AlInSb} & zb & ternary  \\
     \texttt{AlN} & wz & binary \\
     \texttt{AlP} & zb & binary \\
     \texttt{AlPSb} & zb & ternary  \\
     \texttt{AlSb} & zb & binary \\
     \texttt{GaAs} & zb & binary \\
     \texttt{GaAsP} & zb & ternary  \\
     \texttt{GaAsSb} & zb & ternary  \\
     \texttt{GaInP} & zb & ternary  \\
     \texttt{GaInSb} & zb & ternary \\ 
     \texttt{GaN} & wz & binary \\
     \texttt{GaP} & zb & binary 
     \end{tabular}
     \end{minipage}
     \hspace{0.5cm}
     \begin{minipage}[t]{8cm}
     \begin{tabular}{l|c|l}
     \multicolumn{3}{c}{\textbf{Material list}} \\
     \hline
     \textbf{Name}  & \textbf{Crystal structure} & \textbf{Type}  \\
     \hline
     \hline
     \texttt{GaPSb} & zb & ternary \\ 
     \texttt{GaSb} & zb & binary \\
     \texttt{Ge} & zb & simple \\
     \texttt{InAs} & zb & binary \\
     \texttt{InAsP} & zb & ternary  \\
     \texttt{InAsSb} & zb & ternary  \\
     \texttt{InGaAs} & zb & ternary  \\
     \texttt{InGaN} & zb & ternary  \\
     \texttt{InGaSb} & zb & ternary  \\
     \texttt{InN} & wz & binary \\
     \texttt{InP} & zb & binary \\
     \texttt{InPSb} & zb & ternary  \\
     \texttt{InSb} & zb & binary \\
     \texttt{Pentacene} & - & molecular \\
     \texttt{Perowskite} & - & organic  \\
     \texttt{Alq3} & - & organic \\
     \texttt{aNPD} & - & -  \\
     \texttt{Si} & zb & simple \\
     \texttt{SiGe} & zb & binary  \\
     \texttt{SiN} & zb & binary  \\
     \texttt{SiO2} & zb & binary \\
     \texttt{HfO2} & zb & binary \\
     \texttt{TiO2mes} & - & mesoporous \\    
     \texttt{ZnO} & zb & binary \\
     \texttt{ZrO2} & zb & binary 
     \end{tabular}
     \end{minipage}
     \caption{List of Materials}
     \label{table:material_list}
     \end{table}


    
    





The  parameters  of  the  most  important semiconductor materials are  collected  in  the  *material database*.
See :ref:`List of materials <material_list>`      for  a  complete list of  the  materials defined in  tiberCAD. 
*zb* stands  for  *zincblende* crystal  structure, *wz* for  *wurtzite*. 
Materials  can be  simple elements, like Si and  Ge, binary  compounds, such  as GaAs or  GaN,   ternary alloys, like  AlGaAs and quaternary alloys such as AlInGaAs.
For  the  simple and binary  compounds the lattice, strain  and  band properties are included in  each material file. The ternary alloy material  file instead contains  the name of  the  two parent  materials,  from  which the  parameters  of  the  alloy  material are  calculated,  depending  on  the  components concentration, according to a  quadratic  law, e.g. for the energy gap *Eg* 

.. math::
   :label: Vegard_law

    E_g(A_xB_{1-x}) & = xE_g(A)+ (1-x)E_g(B)- x(1-x)C  


where  the  *bowing* parameter *C* accounts  for  the  deviation from a linear interpolation (virtual-crystal approximation) between the two binary compounds *A* and *B* and *x* is the concentration of the binary *A*. Values of *C* for one or  more  parameters are reported in the  alloy  material  file. 
 

The quaternary  alloy material file contains the name of the two parent ternary alloys, e.g. for AlInGaAs ::

  comp_A = GaInAs 
  comp_B = AlGaAs 

along with a default concentration of the two components ::

  x_A = 0.5
  x_B = 0.5

These values can be changed editing the material file.
Instead, when using this quaternary alloy from the input file, one can control the concentration *x* of the *comp_A*, in this case GaInAs.

User-defined  materials can  be  freely  added, provided  that  the  correct syntax for  the parameters entry is  followed.









Input files Examples
----------------------------------


.. _Input_Ex1:

1D Bulk Silicon
^^^^^^^^^^^^



This is the Example 0 in  the  Example directory.
This simple example calculates the IV characteristic of a piece of slightly doped bulk silicon
The mesh  file ``bulk.msh`` has  been  generated by GMSH,  based on the  script bulk.geo,  described in 
:ref:`GMSH_Ex1`. ::

  # Description of the device physical regions
  Device
  {
   meshfile = bulk.msh

   Region bulk 
    {
     material = Si

     Doping
     {
      density = 1e16
      type = donor
     }
    }

  }

Here  the  geometrical and physical structure of the device  to  be  simulated  is  described.
Note that the name of  the  only  ``Region`` of  this  ``Device`` is **bulk**,  that  is  the  name of the Physical Line in  
:ref:`GMSH_Ex1`. ::


  # Definition of Simulation Models 

  Module driftdiffusion
  { 
              
   # name = driftdiffusion  # this is the default name

   #regions = all # 'all' is the default

   # what we want to plot
   plot = (Ec, Ev, eQFermi, hQFermi, ContactCurrent)


   Contact anode { voltage = $Vb }
   Contact cathode { }

  }


Here the block for  **Module**  ``driftdiffusion`` is  defined.

Note  that the names of  the  two ``Contact`` boundary  regions, **anode** and **cathode**, are the  names of the two *Physical Point* objects created in :ref:`GMSH_Ex1`. ::



  # we want to sweep over the anode voltage
  Module sweep
  {
   solve = driftdiffusion
   variable = $Vb
   start = 0.0
   stop = 1
   steps = 10
   # for each step we want to plot the solution variables
   # specified in the driftdiffusion module
   plot_data = true
  }

In  this  example, we  calculate Poisson and  transport equations for a  set  of  bias  values,  as  described  by  **Sweep** block. ::

  Simulation
  {
   # this increases the amount of information
   # written to the screen
   verbose = 2

   solve = sweep

   resultpath = output
   output_format = grace
  }


Finally, the **Simulation**  block defines  the  simulations  to  be  executed and  their  order.
In this  case this  amounts to just the  **sweep**  block calculation.



.. _Input_Ex2:

2D Mosfet 
^^^^^^^^^^^^

This is the Example 4 in  the  Example directory (see  also :ref:`DD_Ex2`)
This  example calculates the IV characteristics of a silicon Mosfet. 
The mesh  file ``mosfet.msh`` has  been  generated by GMSH,  based on the  script *mosfet.geo*,  described in 
:ref:`GMSH_Ex2`. ::

  Device mosfet
  {

   meshfile = mosfet.msh

   material = Si

   Region substrate
   {
     Doping
     {
      density = 1e18
      type = acceptor
     }
   }

   Region contact
   {
     Doping
     {
      density = 5e19
      type = donor
     }
   }

   Region oxide
   {
    material = SiO2
   }

  }


Note that the  ``Regions`` **substrate**, **contact**, **oxide** correspond to the *Physical surfaces* defined in the geo script of GMSH (see  :ref:`GMSH_Ex2`). ::

  Module driftdiffusion
  {

  # we  solve  for electrons only
  #coupling = electrons

  plot = (Ec, Ev, eQFermi, eDensity, eCurrentDensity, eMobility,
          hQFermi, hDensity, hCurrentDensity, hMobility,
          NetRecombination, ElField, ElPotential, ContactCurrents)

  Solver
  {
    type = linesearch

    linear_solver
    {
      method = pconly
      preconditioner = lu
    }
  }


   Physics
   {

    recombination srh {}

    mobility
    {
      type = field_dependent
      low_field_model = doping_dependent
    }
   }

   Contact gate
   {
    type = schottky
    barrier_height = 3.0

    voltage = $Vg

    # assume a gate width of 1 mm = 0.1 cm
    area_factor = 0.1
   }

   Contact source
   {
    voltage = 0.0
    area_factor = 0.1
   }

   Contact backcontact
   {
    voltage = 0.0
    area_factor = 0.1
   }

   Contact drain
   {
    voltage = $Vd
    area_factor = 0.1
   }
  }
             

Note  that the names of  the   ``Contact`` boundary  regions,  are the  names of the  *Physical Line* objects created in :ref:`GMSH_Ex2`. ::





Example of Input file: 2D HEMT
^^^^^^^^^^^^


Here is an example of the input file template::

  Device hemt1
    {
     meshfile = hemt.msh
     mesh_units = 1e-6

     material = GaN
     y-growth-direction = (0,0,0,1)
     z-growth-direction = (1,0,-1,0)

     Region barrier
     { 
       material = AlGaN
       x = 0.14
     }
     Region barrier_doped
     {
       mesh_regions = (barrier_dop_s, barrier_dop_d)
       material = AlGaN
       x = 0.14
       Doping 
       {
         Nd  = 1e21
         type = donor
         doping_level = 0.026
        }
      }

    Region buffer { }
    
     Region buffer_doped
     {
       mesh_regions = (buffer_dop_s, buffer_dop_d)
       Doping 
	{
          Nd  = 1e21
          type = donor
          doping_level = 0.026
         }
      }
     Region cap
     {
       Doping 
       {
         Nd  = 5e18
         type = donor
         doping_level = 0.026
        }
      }
     Region cap_doped
     {
       mesh_regions = (cap_dop_s, cap_dop_d)
       Doping 
      {
        Nd  = 1e21
        type = donor
        doping_level = 0.026
       }
     }
     Region passivation
     {
       mesh_regions = (passiv1, passiv2, passiv3, passiv4)
       material = SiN
     }
     Cluster semiconductor
     {
       regions = -passivation
     }
    }
    Module driftdiffusion
    {
      name = dd
      regions = all
      coupling = electrons
      save_state = true
      #load_state = output_Nt7e16_out/dd_Vd_50.tsv

      plot = (Ec, Ev, Eg, eQFermi, hQFermi, eDensity, hDensity, Polarization,
              eCurrentDensity, hCurrentDensity,CurrentDensity, ContactCurrents,
              ElPotential, ElField, eMobility, hMobility, eConductivity, hConductivity,)

      NonlinearSolver linesearch
      {
        # type = ls
        relative_tolerance = 1e-15
        step_tolerance = 5e-3
        max_iterations = 15
        LinearSolver petsc
        {
          ksp_type = pconly
          preconditioner = lu
        }
      }

      Physics
      {

        band_properties 
          {
             density_of_states bulk_kp 
             {
             strain_simulation = strain
            }
          }

        recombination (srh, direct, auger) { }
        electron_mobility field_dependent
        {
          region = buffer
          low_field_model = constant
         }
        electron_mobility doping_dependent
        {
          region = (barrier, barrier_doped, cap, cap_doped, buffer_doped)
        }
        trap fixed_charge
        { 
          region = GaN/SiN
          Nt = 2.74e13
        }
        trap acceptor
        {
          region = buffer
          Nt = 7e16
          Et = 0.5
          reference = cb
         }
        
        polarization (piezo) {strain_simulation = strain}
        polarization (pyro) {}

      }
	   
     Contact drain
     {
       voltage = @Vd
     }
     Contact source
     {
       voltage = 0.0
     }
     Contact gate
     {
       regions = (gate1, gate2)
       type = schottky
       work_function = 1.5272
       voltage = @Vg
     }
   }
   
 
	
   Module elasticity 
   {

     name = strain
     regions =  -passivation # all the device except Region "passivation"

     plot = (Strain, StrainCell,Stress,Displacement )
     lin_rel_tol = 1e-6 

     Physics 
     {
       body_force lattice_mismatch 
       {
         reference_material = GaN

         structure = wz
         y-growth-direction = (0,0,0,1)
         z-growth-direction = (1,0,-1,0)
          
       }
       
     }

     Contact substrate  
     {type = clamp}

    }


    Module sweep
    {
     name = sweep_g
     variable = Vg
     solve = sweep_d
     start = -1
     stop = 1.0
     steps = 2
     max_step = 0.25
     min_step = 1e-4
     #plot_data = true
    }
	
    Module sweep
    {
     name = sweep_d
     variable = Vd
     solve = driftdiffusion
     start = 0.0
     stop = 50.0
     steps = 200
     min_step = 1e-4
     plot_data = true
    }
	
    Module selfconsistent
    {
     name = relaxation
     solve = (macrostrain, driftdiffusion)
     absolute_tolerance = 1e-3
     relative_tolerance = 1e-5
    }
	
    Simulation
    {
      temperature = 300
      verbose = 3
      solve = (strain, dd, sweep_g)
      logfile = output_Nt5e16_Al0.215_SiN/hemt.log
      # these parameters can be defined in modules, too
      resultpath = output_Nt5e16_Al0.215_SiN
      output_format = vtk
     }
	
	


	
..   </marker>

.. rubric:: Footnotes




