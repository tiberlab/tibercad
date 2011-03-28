..   <marker>

..  _installGetting:


Installation instructions
===================================

In the following, VERSION denotes the version number of the TIBERCAD release you
downloaded and INSTALLPATH denotes the directory where TIBERCAD gets installed.
Version 2.3.0 of **GMSH** (http://www.geuz.org/gmsh) will be installed together with
TIBERCAD. For the Linux version of GMSH you need OpenGL libraries installed on
your system.

Prerequisites
-------------------

Get the installer package for your OS/architecture from http://www.tibercad.org or
by contacting support@tibercad.org. Table 1 lists the packages available for download.
To run TIBERCAD you will also need a license file that you will have to copy into the
installation directory of TIBERCAD.

In the Windows version, some graphical features such as graphical convergence 
monitors are only available if an X Window server is installed and running.

Windows installation procedure
--------------------------------------

To install TIBERCAD in Windows, run the setup program ``tibercad-version_setup.exe`` .

During the installation you can choose the installation directory. 
After finishing installation, copy your license file *(tibercad.lic)* into the 
 **license** subdirectory of the TIBERCAD
installation directory (INSTALLPATH/license), without changing its filename.

| 
| 
| 

..  math::
    :nowrap:
    :label:

    \begin{table}[!h]
    \center
    \begin{tabular}{l|c}
    \multicolumn{2}{c}{\textbf{Installer}} \\
    \hline
    \textbf{installer package name} & \textbf{Target architecture} \\
    \hline
    \hline
     &    \\
    \texttt{tibercad-version_setup.exe } & Windows 32-bit     \\
    \texttt{tibercad-version_installer.bin} & Linux 32-bit self-extracting installer  \\
    \hline
    \end{tabular}
    \caption{Installer Pacakge}
    \end{table}



Linux installation procedure
------------------------------------

To install TIBERCAD in Linux, download and run the self-extracting installer  ``tibercad-version_installer.bin`` 
and follow the installation instructions.

After installation, copy your license file *(tibercad.lic)* into the "license" subdirectory
of the TIBERCAD installation directory (INSTALLPATH/license) without changing the filename. 
You can also provide the license file during installation.

The standard method to launch TIBERCAD is by means of a shell script that is
installed alongside the TiberCAD executable. It takes care of setting all necessary 
environment variables. 

If for some reason you have to run the executable directly, remember
to set TIBERCADROOT to the TiberCAD installation directory (INSTALLPATH).

Quick start guide
----------------------

In the "examples" subdirectory you can find several examples ready to run. They are the
same as the tutorials on http://www.tibercad.org/documentation/tutorial/list.

Windows
^^^^^^^^^^^^^^^^^^^^^

Open Windows Explorer and go to the TIBERCAD installation directory. If you have
write permission in the installation directory, you can browse to an examples directory
and start the simulation by double clicking the input file, e.g. bulk.tib in :ref:`tut0step1` . If
not, copy the whole directory to a location in your personal area and run the examples
from there.

If you cannot run TIBERCAD by double clicking an input file *(\*.tib)* , then the input
files are probably not correctly associated with the TIBERCAD executable. In this case,
try to establish the association by right-clicking the input file, choosing ``open with...`` 

  ``>> Choose Program... >> Browse...`` , browsing to the TIBERCAD installation directory

and choosing the TIBERCAD executable, tibercad.exe. A directory containing the
simulation results will be created with the name provided in the input file.

Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the correct installation of TIBERCAD you should be able to run TIBERCAD from the
command line using the command tibercad. If not, you probably have to add the bin
subdirectory of the TIBERCAD installation directory to your PATH environment variable or
start the TIBERCAD executable using the absolute path (INSTALLPATH/bin/tibercad).

Copy the directory of the example you want to run, e.g. bulk Si, to your home directory
or any place you have write permissions for. Change to the newly created directory and
run TIBERCAD by (assuming :ref:`tut0step1` )

  ``$ tibercad bulk.tib``

A directory containing the simulation results will be created with the name provided in
the input file.

Bug reports / Feedback
-------------------------------------

Please send bug reports, feedback or suggestions to support@tibercad.org. When
submitting bug reports, please always include the full version number of TIBERCAD you
are running. The full version number appears in the first line of output when running
the program:

    $ tibercad

    TiberCAD version 1.0.0-961

    Usage: ``tibercad`` <inputfile>



.. _InputFileGetting:

Input for TiberCAD
=================================================


Input for TIBERCAD is composed by an input file e.g. "input.tib" and a mesh file generated 
by a mesher software: as for now, mesh files from GMSH *(\*.msh, v.1 and v.2.0 )*
and from ISE-TCAD *(\*.grd)* are supported.
Be sure that the material files are in the correct directory .
To run the program, type: **tibercad** *input_file_name*

Description of Input file structure
-------------------------------------------

A valid input file for TIBERCAD is a text file with the structure described in the following.
In the whole input file, everything following a "#" is considered as a comment and is
disregarded; blank lines can be present anywhere and are disregarded too.
Input file is composed by several **blocks** :

A block is enclosed between **"{"** and **"}"** brackets and may include one or more blocks.
Each block has a header made of one or two keywords.
Each block may contain zero or any number of **parameter assignments** in the
form:

     *" tagname = tagvalue"* , where

* *"tagname"* is a string

* *"tagvalue"* is a single numerical or string item or a list of items between "(" and ")"
parenthesis and separated by commas. e.g. ( *cathode, anode* )

Format is free for the parameter assignments, provided that they are separated by
spaces. Everything which follows a "#" is considered as a comment and is disregarded.
For example::

  electron_mobility field_dependent
    {
     region = buffer
     # type = doping_dependent
     low_field_model = constant
    }

Here and in the whole input file a string item can include a combination of characters,
special characters and numbers, except spaces; if a space is found, the string item is
taken as terminated.

The input file is composed by the following main classes of blocks:

    **Device, Module, Simulation**

which will be described in the following.

Device section
-------------------

::

  Device hemt1
    {
     meshfile = hemt_msh.grd
     Region buffer
       {
        .........
        Doping
          {
           .........
          }
       }
     Region barrier_1
       {
        ............
       }
    .................
    }

**Device** section includes the geometrical description of the device to be simulated; an
optional *device name* can be associated to the Device object, e.g. *Device hemt1* .

In **Device** section, two kinds of block can be present: the **Region** block and the
 **Cluster** block. Outside of these blocks, general options common to all the device can
be defined. The most important is the mesh file definition:

  *meshfile* : name of mesh file. N.B.: the extension is mandatory! ( *.grd* for ISE-TCAD,
.msh for GMSH mesh file v.1 and v.2.0 )

The definition of the mesh file associated to this simulation is compulsory:

|  ``meshfile = hemt_msh.grd``
| 

The **Region** blocks contain the description of the device in continuous media approach; 
the **Cluster** blocks define each a group of regions (mesh regions) even with
different physical properties, but to be treated together somewhere in the simulation
(e.g. quantum calculation). In this way it is possible to refer to the set of these regions
simply by the **Cluster** name.

Each **Region** block must be preceded by the keyword **"Region"** , followed by the
(single-word) name of the **TiberCAD Region** . The name of the **TiberCAD Region** 
can coincide with the name of a mesh region, as defined during the modeling of the
device; in this case, if the keyword *mesh_regions* is absent, the **TiberCAD Region** will
be associated to the mesh region identified by the name assigned to the **TiberCAD
Region** . Otherwise, the **TiberCAD Region** will be associated to the mesh regions
specified by the keyword *mesh_regions* .

::

  Region QWell
    {
     mesh_regions = (well1,well2)
     structure = wz
     y-growth-direction = (1,0,-1,0)
     z-growth-direction = (-1,2,-1,0)
     x-growth-direction = (0,0,0,1)
     material = GaN
     Doping
       {
        density = 1e17
        type = donor
        level = 0.025
       }
    }

Here are the description of the available keywords for a **Region** block.

  material (mandatory): name of the material associated to the present region, e.g Si;
it may be a ternary alloy, e.g AlGaAs, in this case keyword x described in the following
has to be present.

  **x** : alloy concentration, expressed as the molar fraction of the first component of
  the alloy; e.g. to express an alloy :math:`Al_xGa_{1-x}As` with molar fraction x = 0.2, that is
  :math:`Al_{0.2}Ga_{0.8}As` , we select **AlGaAs** for the keyword material, and 0.2 for the keyword x,
  thus we write x = 0.2.

  ``mesh_regions`` : (a list of) region physical name(s) as specified in the meshing program.

  ``structure`` : crystal structure (wz = wurtzite, zb = zincblend)

  ``x-growth-direction, y-growth-direction, z-growth-direction`` : Bravais vectors with Miller
indexes for wurtzite crystal (4 element vectors) or zincblende crystal (3 element vectors).

  |warn| A common crystal structure and growth directions definition may be applied to
         all the Regions of the Device, just by defining them everywhere in the Device block, but
         outside any specific Region block.

| 
| 

Subblock doping, with the keywords:

  |  ``density`` : doping concentration [cm-3]
  |  ``type`` : donor or acceptor
  |  ``level`` : energy level of the dopant [eV]

::
	
  Doping
    {
     density = 1e21
     type = donor
     level = 0.026
    }

Each **Cluster** block must be preceded by the keyword **"Cluster"** , followed by the
(single-word) name of the **Cluster** .

::

  Cluster Quantum_1
    {
     regions = (well, barrier1,barrier2)
    }

regions (mandatory): list of physical regions as specified in the meshing program, or
TIBERCAD region names to be grouped in the cluster.

  **Regions** and **Clusters** represent the macroscopical description of the device or structure 
to be be simulated in **TiberCAD** . In the rest of the input file, the physical regions
associated to the **Modules** will be indicated by means of the **TiberCAD Region** and **Cluster** names.

Modules
------------

::

  Module driftdiffusion
    {
     name = dd
     #regions = all
     statistics = FermiDirac
     strain_simulation = macrostrain
     plot = (Ec, Ev, Eg, eQFermi, hQFermi, eDensity, hDensity, Polarization,
     eCurrentDensity, hCurrentDensity,CurrentDensity, ContactCurrents,
     ElPotential, ElField, eMobility, hMobility, eConductivity, hConductivity)
      ..........
     Solver
       {
        nonlinear_solver = tiber
        nonlin_max_it = 15
        nonlin_rel_tol = 1e-12
        #nonlin_abs_tol = 1e-15
        nonlin_step_tol = 1e-5
       }
     Physics
       {
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
the keyword **"Module"** , followed by the (single-word) module name. This must be the
name of one of the TIBERCAD modules.
Here are the Modules implemented until now:

| ``driftdiffusion`` : Poisson-driftdiffusion transport of electrons and holes

| ``thermal`` : Heat balance simulation

| ``excitontransport`` : Exciton transport model

| ``macrostrain`` : Calculation of Elastic deformations in heterostructures

| ``efaschroedinger`` : Envelop Function Approximation (EFA) solution of single particle 
Schroedinger equation for electrons and holes

| ``quantumdensity`` : Calculation of quantum density of electrons and holes.

| ``quantumdispersion`` : Dispersion of quantized states in k space

| ``opticskp`` : Optical properties (optical kp matrix elements)

| ``opticalspectrum`` : Emission spectrum (with k-space integration)

| ``Sweep`` : Parameterized execution of a module simulation (e.g. for the calculation
of output current characteristics)

| ``Selfconsistent`` : coupled calculations of different simulation modules.

| ``DSC`` : Simulation of a DSC solar cell1.

Each module-block usually contains a list of general options, such as plot and others
specific to each module. Then, two main blocks define the Physics and the Solver
models and parameters for this module.

  **Solver** contains the solver parameters; depending on the Module, it can contain a
LinearSolver definition subblock.

  **Physics** usually contains the definition of the physical models used in the simulation.


For example, for driftdiffusion module, recombination, electron mobility, trap, polarization 
and so on. A particular model is the Boundary model, which has an alias Contact
for driftdiffusion module.
The declaration of these models obey to the following syntax:

    *model_keyword* type_specifier
    <block>

where *model_keyword* is the name of the physical model to be declared (e.g. trap
model) and **type_specifier** is the name of a particular one among the available descriptions 
for that model. For example, for a trap model of **type** *acceptor* ::

  trap acceptor
    {
     region = buffer
     Nt = 7e16
     Et = 0.5
     reference = cb
    }

An alternative multiple declaration is possible if no parameters, beyond default, are
declared:

    *model_keyword* (type specifier1, type specifier2,...)

  ``recombination (srh, direct, auger) { }``

In this example, several recombination models are defined (srh, auger, direct ) each one
with default parameters.
For a detailed description of the models, please refer to the reference guide.
Two special Modules are: the Sweep Module and the Selfconsistent Module

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
  Module sweep
    {
     name = sweep_gate
     solve = sweep_drain
     variable = $Vg
     start = -0.1
     stop = 1.0
     steps = 11
    }

Each sweep Module defines a set of calculations applied to a boundary region (e.g. a
set of bias values to be assigned to a drain contact of a MOSFET for the calculation of
an output drain IV characteristic), in this Guide referred to as sweep calculation.

The following keywords are defined for this feature:

|  ``variable`` : name of the variable to which the sweep is applied: 

E.g.: variable = $Vg

indicates that the values are applied to the variable Vg, which is a quantity defined in a **Contact** definition

|  ``start, stop, steps`` : sweep starts from start value, is repeated steps times and stops
in stop

|  ``solve`` : name of the simulation (module) associated to the sweep calculation; it may
be the name of another sweep defined in the same block.

|  ``plotvariable`` (obsolete): specify the integrated quantity to be calculated during the
sweep and that will be shown in the output file ``sweep_modelname_sweepvariable.dat`` ,
eg. sweep driftdiffusion Vb.dat for a sweep of current calculation on the variable Vb
(typically a contact voltage).

|  ``plot_data`` : default is false; if it is set to true, then output data will be written for
each step of the sweep calculation, otherwise just the results for the final step will be
present in the output.

Once a *sweep* calculation has been defined, it is treated as a special case of simulation
and may be executed as an usual simulation: by adding it in the solve list, 
e.g. ``solve = sweep_drain`` .

Module selfconsistent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this Module it is possible to define a self-consistent calculation based on two different
simulation modules (e.g. driftdiffusion and excitontransport).

::

  Module selfconsistent
    {
     name = sc_all
     solve = (tb, dens_el, dens_hl, dd)
     max_iterations = 5
     abs_tolerance = 1e-3
     rel_tolerance = 1e-6
    }

In **solve** the list of simulations to be performed self-consistently is specified.
Now it is possible to execute the specified simulations in self consistent way, by using
the **selfconsistent** keyword like a simulation name, in any **solve** assignment, e.g. 
``solve = selfconsistent`` in **Simulation** section, or even in a **sweep** section.

Simulation section
--------------------------

In this block one can specify several general parameters and settings for the actual
calculation to be run, such as the temperature, the process-flow of simulation, etc.

|  ``searchpath`` : path for material files

|  ``mesh units`` : units of measurements used in the meshing (relative to meters): e.g.,
:math:`10^{-6}` for :math:`\mu m`

|  ``dimension`` : dimension of simulation (1,2,3)

|  ``temperature`` : temperature of the system [K]

|  ``solve`` : list of simulations to be executed, in the order of execution; if the list contains
"sweep", a sweep is performed as specified in sweep block in the Solver section.

|  ``solve = (strain,driftdiffusion, quantum_electrons, quantum_holes)``

|  ``resultpath`` : path for output directory

|  ``output format`` : format of the output data: *gmv* for **GMV** , ise for **Tecplot** , grace for
xmgr (ascii data column type), *vtk* for **Paraview** .

|  ``plot`` : list of output variables which are calculated and available in output files. It
may be overriden by defining list specific to one or more modules.

Output description
-------------------------

At the end of the execution, the program will write the results of the simulation in the
directory specified by resultpath , with the format specified by output format. The output
variables are specified in a list plot, in each Module.

TiberCAD output is divided in two classes: **mesh-based** and **mesh-independent**
quantities.

**Mesh-based** quantities are all the quantities associated with the nodes of the mesh,
such as Fermi level, electron and hole density, conduction and valence band, etc., together
with all the quantities associated with the elements of the mesh, such as current density.

The output values for these quantities are reported in the files simname msh.ext,
where *simname* is the simulation module used for the calculations and ext is the extension
of the chosen file format, e.g. vtu for paraview output.

|  ``strain_msh.vtu``

In the case a sweep calculation is performed and the plot data keyword is set to
true, the output files are of the kind simname sweepvariable step value msh.ext, where
 ``sweepvariable`` is the variable with respect to which the sweep is performed (e.g. gate
voltage) and step value is the value of this variable at that step; e.g the result at the
step *Vbias = 1.1* will be found in the file:

|  ``dd_Vbias_1.1_msh.vtu``

|  ``mesh-independent`` quantities are the quantities which are not associated to the
mesh, for example current at the contacts of a diode or quantized energy levels in a
quantum well. These **mesh-independent** quantities are displayed in separated files,
with the format simname.ext, e.g ``quantum_electrons.dat`` , where simname is the name of
the model (simulation) associated to the results. If a sweep is performed, the output file
gets the format sweep name simname.ext, where *sweep_name* is the name of the sweep
performed, for example

|    ``sweep_drain_driftdiffusion.dat``

Inside the file, output values for all the steps of calculation are shown.

Example of Input file
----------------------------------

Here is an example of the input file template::

  Device hemt1
    {
     meshfile = hemt_msh.grd
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
           density = 1e21
           type = donor
           level = 0.026
          }
       }
    Region buffer { }
    
	Region buffer_doped
      {
       mesh_regions = (buffer_dop_s, buffer_dop_d)
       Doping 
	     {
          density = 1e21
          type = donor
          level = 0.026
         }
      }
    Region cap
      {
       Doping 
	     {
          density = 5e18
          type = donor
          level = 0.026
         }
      }
    Region cap_doped
      {
       mesh_regions = (cap_dop_s, cap_dop_d)
       Doping 
	     {
          density = 1e21
          type = donor
          level = 0.026
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
   
  Module macrostrain
    {
     regions = -passivation # all the device except Region "passivation"
     plot = all
     LinearSolver # default
       {
        tolerance = 1e-7
        max_iterations = 10000
        #xmonitor = true
        pc = composite
       }
     Physics
       {
        Boundary substrate
          {
           regions = bottom
           material = GaN
           y-growth-direction = (0,0,0,1)
           z-growth-direction = (1,0,-1,0)
           x-growth-direction = (-1,2,-1,0)
          }
		  
        Boundary extended_material
          {
           regions = side
          }
       }
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
     solve = (macrostrain, dd, sweep_g)
     logfile = output_Nt5e16_Al0.215_SiN/hemt.log
     # these parameters can be defined in modules, too
     resultpath = output_Nt5e16_Al0.215_SiN
     output_format = vtk
    }
	
	
.. rubric:: Footnotes

	
..   </marker>


.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%

