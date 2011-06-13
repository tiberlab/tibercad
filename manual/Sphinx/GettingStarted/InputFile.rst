..   <marker>

..  _installGetting:


Installation instructions
===================================

In the following, ``installpath`` denotes the directory where tiberCAD 2.0.0 gets installed.
Version 2.5.0 of ``GMSH`` (http://www.geuz.org/gmsh) will be installed together with
tiberCAD.

Prerequisites
-------------------

Get the installer package for your OS/architecture from http://www.tiberlab.com or
by contacting support@tiberlab.com. Table :ref:`Installer packages<input_installerpacks>`  lists the packages available for download.
To run tiberCAD you will also need a license file that you will have to copy into the
installation directory of tiberCAD.

.. In the Windows version, some graphical features such as graphical convergence 
.. monitors are only available if an X Window server is installed and running.

Windows installation procedure
--------------------------------------

To install tiberCAD in Windows, run the setup program ``tibercad-2.0.0_setup.exe`` .

During the installation you can choose the installation directory. 
After finishing installation, copy your license file ``tibercad.lic`` into the ``license`` subdirectory of the tiberCAD
installation directory (``installpath/license``), without changing the filename.


.. _input_installerpacks :

..  math::
    :nowrap:

    \begin{table}[!h]
    \center
    \begin{tabular}{l|c}
    \multicolumn{2}{c}{\textbf{Installer}} \\
    \hline
    \textbf{installer package name} & \textbf{Target architecture} \\
    \hline
    \hline
     &    \\
    \texttt{tibercad-2.0.0-setup.exe} & Windows 32-bit     \\
    \texttt{tibercad-2.0.0-installer.bin} & Linux 32-bit self-extracting installer  \\
    \hline
    \end{tabular}
    \caption*{Installer Packages}
    \end{table}



Linux installation procedure
------------------------------------

To install tiberCAD under Linux, download and run the self-extracting installer
``tibercad-2.0.0-ARCH_installer.bin``, where ARCH corresponds to your hardware architecture,
and follow the installation instructions.

After installation, copy your license file ``tibercad.lic`` into the ``license`` subdirectory
of the tiberCAD installation path (``installpath/license``) without changing the filename. 
You can also provide the license file during installation.

tiberCAD is launched by means of a shell script, which is
installed together with the tiberCAD executable. It takes care of setting all necessary 
environment variables. 

If you have to call the executable directly, you have to
set the environment variable ``TIBERCADROOT`` to the tiberCAD installation directory (``installpath``).

Quick start guide
----------------------

In the ``examples`` subdirectory you can find several examples ready to run. More examples may be
available on http://www.tiberlab.com or http://www.tibercad.org.

tiberCAD has the following command line options:

  -v     print the version number and exit
  -b     run in batch mode, without asking for user input.

The ``-b`` option is useful only for the windows version, which attends a keystroke by the user before
exiting.

Windows
^^^^^^^^^^^^^^^^^^^^^

Open Windows Explorer and go to the tiberCAD installation directory. If you have
write permission in the installation directory, you can browse to an example directory
and start the simulation by double clicking the input file, e.g. ``bulk.tib`` in :ref:`tut0step1` . If
not, copy the whole directory to a location in your personal area and run the examples
from there.

If you cannot run tiberCAD by double clicking a tiberCAD input file (the files with extension ``.tib``), then the input
files may not be correctly associated with the tiberCAD executable. In this case,
try to establish the association by right-clicking the input file, choosing 

  ``open with... >> Choose Program... >> Browse...``,
browsing to the tiberCAD installation directory and choosing the tiberCAD executable, ``tibercad.exe``.

Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After a correct installation of tiberCAD you should be able to run tiberCAD from the
command line using the command ``tibercad``. If not, you probably have to add the ``bin``
subdirectory of the tiberCAD installation path to your ``PATH`` environment variable or
start the tiberCAD executable using the absolute path (``installpath/bin/tibercad``).

Copy the directory of the example you want to run to your home directory
or any place you have write permissions for. Change to the example directory and
run tiberCAD by (assuming :ref:`tut0step1` )

  ``$ tibercad bulk.tib``


Bug reports / Feedback
-------------------------------------

Please send bug reports, feedback or suggestions to support@tibercad.org. When
submitting bug reports, please always include the full version number of tiberCAD you
are running. The full version number appears in the first line of output when running
the program:

::

  $ tibercad -v
  tiberCAD version 2.0.0 (x86_64-linux)



.. _InputFileGetting:

Input for tiberCAD
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
In the ``Device`` section, two types of blocks can be present: ``Region`` and 
``Cluster`` blocks. Outside of these blocks, general options common to all the device can
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
       Nd  = 1e17
      type = donor
      doping_level  = 0.025
    }
  }

In  the following, instead, the  ``Region``  *well1* takes  its  name  from the mesh  region *well1*  it  refers to

::

  Region well1
  {
    
    Doping
    {
       Nd  = 1e17
      type = donor
      doping_level  = 0.025
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
          Bravais vectors with Miller indices for wurtzite (4-tuple) or zincblende (3-tuple) crystal
          along the x, y and z directions.

.. ``structure`` : string
..          crystal structure (wz = wurtzite, zb = zincblend)


.. tip:: A common material, crystal structure and growth directions can be defined for all
          device regions by defining them outside of the ``Region`` blocks.



The optional subblock doping as in the example above contains the keywords:

  ``type`` : string
         The dopant type. Can be ``donor`` or ``acceptor``.

  ``Nd`` : double 
         The doping concentration in cm\ :sup:`-3`.

  ``doping_level`` : double
         The energy level of the dopant given as the distance from the conduction band edge (for donors)
         or from the valence band edge (for acceptors) in eV.

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
      strain_simulation = strain

      polarization (piezo, pyro) {}
      
      mobility
      {
        type = field_dependent
        low_field_model = doping_dependent
      }
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

 ``quantumdispersion`` : 
       Dispersion of quantized states in k space

 ``opticskp`` : 
       Optical properties (spontaneous emission  spectrum in k=0)

 ``opticalspectrum`` : 
       Emission spectrum (with k-space integration)

 ``DSC`` : 
       Simulation of a DSC solar cell

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

Writing e.g. ``variable = $Vg`` that a ``sweep`` on the variable ``$Vg`` will be performed.

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
     If set to true, output data will be written after
     each step of the ``sweep``, otherwise this will be done only for the final sweep value.
     The default is ``false`` to limit output data.
     The amount of written data can be controlled using the ``max_step`` option

Once a ``sweep`` calculation has been defined, it is treated as an ordinary simulation
and may therefore be used like any simulation by adding it to any ``solve`` option.

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

  ``output_format`` : 
     format of the output data: ``grace`` for ascii data column type
     ``vtk`` for VTK based output.


Output description
-------------------------

At the end of the execution, the program will write the results of the simulation in the
directory specified by  ``resultpath``, with the format specified by ``output_format``. 
The output variables for each Module are specified in the list *plot*.

tiberCAD output is divided in two classes: **mesh-based** and **mesh-independent**
quantities.

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
gets the format ``sweepname_simname.ext``,  e.g. ``sweep_drain_driftdiffusion.dat``, where ``sweep_name`` is the name of the ``sweep``
performed, for example

  ``sweep_drain_driftdiffusion.dat``
The mesh independent results for every sweep step are stored in this file.



Example of Input file
----------------------------------

Here is an example of the input file template::

  Device hemt1
    {
     meshfile = hemt.msh
     mesh_units = 1e-6

     material = GaN
     y-growth-direction = (0,0,0,1)
     z-growth-direction = (1,0,-1,0)
     x-growth-direction = (-1,2,-1,0)

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
      statistics = FermiDirac
      strain_simulation = macrostrain
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
        polarization (piezo, pyro) { }
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
         x-growth-direction = (-1,2,-1,0)
          
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




