.. highligthlang:: c

.. _overview:

Overview
========

Introduction to numerical simulation with Tiber-CAD
--------------------------------------------------


TiberCAD is a multiphysics software tool; it comprises a set of solvers, called simulation
models, each one describing a physical problem to be solved, e.g. DriftDiffusion (to solve
Poisson and DriftDiffusion equations), EFASchroedinger (to solve Schroedinger equation
in envelope function approximation), Macrostrain (to calculate macroscopical strain with
an elastic model) and others.

Similarly to other device CADs, TiberCAD requires that one follows a three-step
procedure. [#]_

In a first step the device geometry must be sketched, giving all the geometrical information 
needed by the simulations. This can be performed by means of a text  file
or with the help of a graphical tool. During this procedure one or more *mesh_regions*
and *boundary regions* have to be defined: in a following stage the *mesh_regions* will be
associated to materials and *boundary regions* to device contacts (boundary conditions in
general).

The second step consists in running a *mesher* tool, which reads the geometry  file
and sets up the computational mesh used to discretize the partial differential equations
representing the physical models to be solved. For this procedure TIBERCAD makes use
of the GPL software *gmsh*. Optionally, mesh output provided by other meshers, such
as the one included in ISE-TCAD tool, are also supported. The output of the meshing
procedure is a mesh  file that contains information about the space discretization as well
as the *mesh_regions* and the *boundary regions*.

In the last step the actual simulations are performed. Together with the mesh information 
(comprised in the mesh  file ), TiberCAD requires an input  file which associates 
*materials* to *mesh_regions*, defines the physical properties and physical models to be 
applied, and the type of calculations to be performed.

|more| Details about the modeler and mesher tools can be found in the specific user-manuals.
Here we deal primarily with the TiberCAD input  file. However, in discussing examples
of 1D, 2D and 3D simulations, we will also describe in some detail the geometry input
files used to run *gmsh*.
 

Input file structure
--------------------

TiberCAD input  file is a text  file which includes a description of the device structure,
the definition of the solvers to be executed, with all the relevant physical and numerical
parameters for each of them.

The input  file is organized into several sections, each describing a different aspect of
the problem to be solved. The strategy employed is similar to other commercial T-CAD
tools and requires some practice to reach a good level of familiarity. We strongly suggest
to read first the following chapters ( :ref:`GettingStarted1D` and :ref:`GettingStarted2D` ) in this manual and then
study the input  files provided in the examples directory and try to modify them, before
writing your own input  file from scratch. The example  files touch all current features
implemented in the code.

Let's see an overview of the main features of the input  file The core of the input  file
comprises three sections, called **Device**, **Models** and **Simulation**::

  $Device
    {
     Region Si_channel
       {
        material= Si
        doping = 1e16
       }
     Region gate_oxide
      {
       material= SiO2
      }
    ...
    }
  
  $Models
    {
     model driftdiffusion
      {
       simulation_name= dd
      }
    ...
    }
  $Simulation
    {
     solve= dd
    }
  
The section **Device** is used to associate one or more mesh_regions to a material and
to a set of physical properties such as doping concentrations, doping levels, etc. The
section **Model** is used to define the physical models to be solved. Each physical model
can be applied to the whole device, or to a set of regions of the device, (defined in **Device**
section). Finally, **Simulation** section states the type of calculation to be executed, that
is the *simulation* to be *solved*. These sections will be described in full details in the
following.

Besides these main sections, there other other two, called **Solver** and **Physics**, where
some parameters can be set, respectively for the numerical solvers and for the physical
models. The aim of these sections is to give the user the maximum of flexibility 
to tune his/her simulation. Especially regarding the Solver case, values of the numerical
parameters have been already tuned for each application and should be modified only
by an advanced user.

.. _DefinitionRegions: Boundary

Definition of physical and boundary regions in TiberCAD
-------------------------------------------------------

Let's see now in more detail how to associate *physical* information to the regions of the
device model and how to define the boundary regions, that is contacts or, in general,
regions where some kind of boundary condition is to be applied.

When TIBERCAD is run, it reads the mesh file which contains the finite element grid
which meshes the geometrical description of the device or nanostructure, and which will
be the basis of PDE discretization.

As we have seen before, to execute the proper simulations, TIBERCAD needs some in
formation about the physical and boundary regions associated with the mesh. A physical 
region associates all the elements corresponding to an homogeneous part of the device
(usually related to the same material or doping). In TIBERCAD, these regions are referred
to as **mesh_regions**.

As for boundary regions, they are needed to specify boundary conditions (b.c.) for
the solution of the PDEs of our simulation. By default, to all the external boundary of
the device a Neumann b.c. is imposed, meaning null derivative of electric field and zero flux 
of current normal to the boundary. These are the usual b.c. applied in the simulation
of electronic devices; in particular, these conditions are implicitly satisfied by using the
finite element formulation.. Usually, however, one needs to impose also specifical b.c.
to the device, relative, most often, to contacts of some kind (ohmic, schottky), but also
heat and temperature b.c. or reference substrates for strain calculations. These regions,
(constituted by surfaces, lines or points, respectively for 3D, 2D and 1D simulations) are
called in TIBERCAD **boundary regions**.

It is important to know that the information about the physical and boundary regions
must be present in the mesh file before it is read by TIBERCAD, and thus have to be
produced by making use of the modeling/mesher software. As for now, TIBERCAD
supports the mesh output of the following software tools: **GMSH v.2** and **ISE-TCAD
v.9.5**.

Using ISE-TCAD
^^^^^^^^^^^^^^

By means of the utility DEVISE of ISE-TCAD v.9.5, it is possible to design and mesh
a device; after the meshing has been successfully performed, an output file is produced,
with the extension *.grd*. This file contains the description of the mesh and also the
list of the user defined material regions and contact regions. By reading this *.grd* file
in TIBERCAD, one can refer to the ISE TCAD material regions, simply with the user
defined name, which is present in the *.grd* file. This name should be unique in the whole
device. In the same way, ISE TCAD Boundary regions (*Contacts*) can be referred to in
TIBERCAD by means of their user-defined name, present in the *.grd* output file, too.

Using GMSH
^^^^^^^^^^

If GMSH program is used to model and mesh the device, a bit more care has to be
taken. Here we introduce the procedure to be followed; in the following tutorials ( :ref:`GettingStarted1D` )
The subject will be considered in detail with a step-by-step description. See
also the GMSH user manual (http://geuz.org/gmsh) for further details.

In the context of GMSH, it is possible to define several 1, 2 and 3D *Physical Entities*.
These *Physical Entities* allow to associate one or more geometrical entities to a single
numerical ID or, better, to a string label, so that several **mesh_regions** and **boundary
regions** can be defined and referred to by TiberCAD. Here is a simple example of a
script to generate a 1D geometrical model ( *\*.geo*) file in GMSH ( see also :ref:`GettingStarted1D` )::

  Point(1) = {-25,0,0,0.5};
  Point(2) = {0,0,0,0.002};
  Point(3) = {25,0,0,0.5};
  Line(1) = {1,2};
  Line(2) = {2,3};

  Physical Line("p_side") = {2};
  Physical Line("n_side") = {1};
  Physical Point("cathode") = {3};
  Physical Point("anode") = {1};

Here, first the geometrical entities **Points** and **Lines** are defined.

In the definition of the geometrical **Points**, the three first expressions inside the
braces on the right hand side give the three X, Y and Z coordinates of the point ; the
last expression (0.5 or 0.002 in this example) sets the **characteristic mesh length** at
that point, that is the **size** of a mesh element, defined as the length of the segment for
a line mesh element, the radius of the circumscribed circle for a triangle mesh element
and the radius of the circumscribed sphere for a tetrahedron mesh element. Thus, the
smaller is the value of the **characteristic mesh length**, the greater is the mesh density
close to that point. The size of the mesh elements will then be computed in GMSH by
linearly interpolating these characteristic lengths in the whole mesh.

In the definition of a geometrical **Line**, the two expressions inside the braces on the
right hand side give the identification numbers of the start and end **Points** of the line.
Then, two physical regions are defined, each associated to one of the two geometrical
entities: **Physical Line(n side)** and **Physical Line(p side)**. The expression(s) inside
the braces on the right hand side give, in general, the identification numbers of all the
geometrical lines that need to be grouped inside the **Physical Line**.
In this way, these physical regions are made available for TIBERCAD, and will be used
to associate them to a TIBERCAD region through the keyword **Region**, as follows::

  Region n_side
    {
      ............
  ............
  Region p_side
    {
     .............
  ...........

It is also possible to group more than one physical region in a single Device Region,
with the keyword mesh_regions, as follows::

  Region reg_1
    {
     mesh_regions = (region1, region2)
     ............
  Region reg_2
    {
     mesh_regions = (region3, region4)
  ...........

Note that, in this case, *region1, region2, region3, region4* are the labels previously
defined inside GMSH, while *reg_1* and *reg_2* are names chosen by user for his convenience.

Then, in the GMSH script, two **Physical Point** are defined, **anode** and **cathode**,
and associated to the first and to the last point of our 1D device. These points are
needed to impose some boundary conditions and in this way they are made available for
TiberCAD, and will be used to associate each of them to a boundary condition region,
through the keyword **BC_Region**::

  BC_Regions
    {
     BC_Region cathode
       {
        .........
        ......................
     BC_Region anode
       {
        ........
        ....................
  
As before, it is possible group more boundary regions in a single **BC_Region**, with
the keyword **BC_reg_numb**::

  BC_Regions
    {
     BC_Region cathode
       {
        BC_reg_numb = (contact1, contact2)
        ......................
     ....................

Again, in this second case, *cathode* is a label chosen by user, while *contact1* and
*contact2* are labels previously defined inside GMSH.

In 2D case, a set of **Physical Surface** will be defined to be used as **mesh_regions**,
while **Physical Line** is used for **boundary regions**.

Finally, in 3D case, **Physical Volume** is used to define **mesh_regions**, while **Physical 
Surface** is used to define **boundary regions**.

Simulation environments
-----------------------

TIBERCAD allows to compute different physical models in different parts of a device or
nanostructure by coupling in a general way different *simulation environments*. A *simulation 
environment* is composed by all the physical regions to which a particular model is
assigned. A *simulation environment* is therefore defined by the mesh elements belonging
to its physical regions and by the **simulation model** which has been associated to these
regions. This association is made possible by the definition of TIBERCAD **Regions** and
**Clusters**. Different **simulation environments** can have physical regions in common.

In this way, each simulation is run on a subset of the device and can be possibly coupled
(even self-consistently) with a simulation run on a different subset of the device corresponding 
to a *simulation environment* with a non-void intersection with the first one.

The coupling of the two simulations is performed by means of appropriate Boundary
Conditions (e.g. Current Density, Voltage, ...). In principle, the two *simulation environments* 
can refer to two simulations at different scale, e.g. atomistic tight-binding and
macroscopic drift-diffusion. This allows an effective **multi-scale** simulation of the device
to be studied.


.. |more| image:: more.*
          :scale: 50%

.. rubric:: Footnotes

.. [#] example of Footnotes

