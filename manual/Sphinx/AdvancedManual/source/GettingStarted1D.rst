.. highligthlang:: c

.. _GettingStarted1D:

Getting started 1D
==================


In this section we will see, step by step, how to use TIBERCAD to simulate numerically
a semiconductor device. 

|more| As a very simple example we will refer to the Tutorial_00_ (Si
bulk) that you can find in the :ref:`Tutorial` chapter.

.. _Step1D: Steps

..  index:: pair:Device;modeling

Step 1: Modeling the device
---------------------------

As a first step, we have to model the device. To do so, you can use DEVISE module
of ISE-TCAD 9.5 software package or GMSH program. Here we'll see in details the
procedure for GMSH. There are two possible ways to use GMSH:

  #. Interactive, using the graphical interface
  
  #. Using a script file.

In the following we'll see how to write a basic GMSH script (bulk.geo); for any details
please refer to GMSH manual GMSH (http://geuz.org/gmsh/).
In a GMSH script, several variables can be defined and given a value in this way::

  L = 1;
  d = 0.01;

these are valid GMSH variables: *L* is just the length of the Si sample; *d* is the value
of a **characteristic mesh length** (see below).

     * Definition of geometrical entities Points::

  Point(1) = {0, 0, 0, d};
  Point(2) = {L, 0, 0, d};

In the definition of a geometrical point, the three first expressions inside the braces
on the right hand side give the three X, Y and Z coordinates of the point; the last
expression *(d)* sets the **characteristic mesh length** at that point, that is the size of a
mesh element, defined as the length of the segment for a line mesh element, the radius of
the circumscribed circle for a triangle mesh element and the radius of the circumscribed
sphere for a tetrahedron mesh element.

Thus, the smaller is the value of *d*, the greater is the mesh density close to that point.
The size of the mesh elements will then be computed in GMSH by linearly interpolating
these characteristic lengths in the whole mesh.

|warn| In a 1D simulation it is assumed that the geometrical model is restricted to
the x axis. Any other geometrical orientation could give unpredictable results.

  * Definition of geometrical entity **Line**::

  Line(1) = {1, 2};

The two expressions inside the braces on the right hand side give the identification
numbers of the start and end points of the line.


  * Definition of the physical entity **Physical Line 1**::

  Physical Line(1) = {1};

The expression(s) inside the braces on the right hand side give the identification
numbers of all the geometrical lines that need to be grouped inside the **physical line**.
In this way, in general, **physical regions** are created which associate together geometrical 
regions, and then the related mesh elements, which share some common physical
properties. It's only these physical regions which can be referred to outside GMSH. In
TIBERCAD, this is done by associating one or more physical regions to a **TiberCAD**
region through the keyword **mesh_regions** (see in the following).

  * Definition of two physical entities **Physical Point**::

  Physical Point(1) = {1};
  Physical Point(2) = {2};

Beginning from GMSH v.2, it is possible, alternatively, to assign more convenient
**Physical Names** to the Physical entities, instead of numerical IDs. **Physical Names**
consist of strings enclosed between quotation marks. The syntax is the following::

  Physical Line("bulk") = {1}
  ............
  Physical Point("Anode") = {1};

.. note::

  |idea| In general, in a n-Dimension *(nD)* simulation, *(n-1)D* physical regions (points
  in 1D, lines in 2D, surfaces in 3D) are used by TIBERCAD to impose the required boundary
  conditions. Each *(n-1)D* physical region defined in this way in GMSH will be associated
  in TIBERCAD to a boundary condition region, through the keyword **BC_reg_numb**.
  Thus, in this case, Physical points **1** and **2** will be associated respectively to two BC
  regions (see in the following).

..  index:: pair:Device;meshing
  
Step 2: Meshing the device
--------------------------

The *.geo* script file with the geometrical description can be run in GMSH, to display the
modelled device and to mesh it through the GMSH graphical interface. Alternatively, a
*non-interactive* mode is also available in GMSH, without graphical user interface. For
example, to mesh this 1D tutorial in non-interactive mode, just type::

    gmsh bulk.geo -1 -o bulk.msh
  
where bulk.geo is the geometrical description of the device with GMSH syntax; 

    -1 means 1D mesh generation;

some command line options are::

  -1, -2, -3 to perform 1D, 2D or 3D mesh generation,
  -o mesh file.msh to specify the name of the mesh file to be generated

In this way, a *.msh* has been generated and is ready to be read in TIBERCAD.

..  index:: pair:Device;input file

Step 3: TiberCAD Input file
------------------------

Now we have to write down the TIBERCAD **input file** ( see *bulk.tib* in the :ref:`Tutorial` ).

Definition of Device Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we have to list all the TIBERCAD Regions present in our device: a TIBERCAD
**Region** is usually a somehow physically homogeneous region of the device or the nanostructure 
we are going to model, featuring the same material and possibly the same
doping. ::

  $Device
    {
     Region bulk
       {
        mesh_regions = 1
        material = Si
        doping = 1e16 doping_type = donor
      }
    }
    
In this example, the TIBERCAD **Region** bulk is made of Silicon and n-doped with a
concentration :math:`10^{16}cm^{-3}`.

Through the keyword **mesh_regions**, one or more of the *physical regions* (*Physical
Lines* in 1D, *Physical Surfaces* in 2D, *Physical Volumes* in 3D) previously defined in the
GMSH mesh can be associated to the present TIBERCAD **Region**.

With *mesh_regions = 1*, we associate the Physical Line 1, defined in the Step 1, to
the TIBERCAD **Region** bulk.

Definition of Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we define the Simulation *driftdiffusion_1* : it belongs to the class **driftdiffusion**::

  Models
    {
     model driftdiffusion
       {
        options
          {
           simulation_name = driftdiffusion_1
           physical_regions = all
          }

The TIBERCAD simulation *driftdiffusion_1*, belonging to the **model driftdiffusion**,
will be applied to the whole device structure (*physical_regions = all*).

Definition of Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The anode and cathode contacts of our 1D Si sample are defined as **Boundary conditions 
regions** (*BC_Region anode, BC_Region cathode*) in the following way::
  
  BC_Region anode
    {
     BC_reg_numb = 1
     type = ohmic
     voltage = @Vb
    }
  BC_Region cathode
    {
     BC_reg_numb = 2
     type = ohmic
     voltage = 0.0
    }
 
Both contacts are defined as *ohmic*, cathode is assigned a fixed *voltage = 0.0*, while anode
voltage is given by the value of the variable *Vb (voltage = @Vb)*.

Through the keyword **BC_reg_numb**, one or more of the **(n-1)-Dimension** physical
regions (Physical Points in 1D, Physical Lines in 2D, Physical Surfaces in 3D) previously
defined in the GMSH mesh can be associated to the present TIBERCAD **BC Region**.

With *BC_reg_numb = 1*, we associate the Physical Point 1, defined in the **Step 1**, to the
TIBERCAD **BC_Region** *anode*; with *BC_reg_numb = 2*, we associate the Physical Point
2, defined in the **Step 1**, to the TIBERCAD **BC_Region** *cathode*.

Alternatively, one can make use of the **physical names** associated to the physical
regions in the meshing tool. In this case, we simply associate the **(n-1)D** physical region,
respectively *"anode"* and *"cathode"* , by means of the TIBERCAD **BC_Region name**::

  BC_Region anode
    {
     type = ohmic
     voltage = @Vb
    }
  BC_Region cathode
    {
     type = ohmic
     voltage = 0.0
    }
  
Note that, in this case, the TIBERCAD **BC Region name** needs to be identical to
one of the **physical names** defined during the modeling of the device with GMSH.

Definition of Simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The variable *Vb* is specified in the *sweep* block, in the **Solver** section::

  sweep
    {
     simulation = driftdiffusion_1
     variable = Vb
     start = 0.0
     stop = 1
     steps = 10
    }
  
In this way, the simulation *driftdiffusion_1* is performed for 10 (*steps = 10*) values
of the anode voltage (*variable = Vb*), between 0 and 1.

Definition of Execution parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the **Simulation** section, we decide *which* simulations to perform and in which *order*;
 we set *solve = sweep*, to execute the sweep which run *driftdiffusion_1* **simulation** for
the specified loop::

  $Simulation
    {
     # searchpath = .
     meshfile = bulk.msh
     dimension = 1
     temperature = 300
     solve = sweep
     resultpath = output
     output_format = grace
     plot = (Ec, Ev, ContactCurrents)
    }
  
Output files with conduction and valence band profiles *(plot = Ec,Ev..)* and all the
calculated values of the current at the contacts *(ContactCurrents)* (the IV characteristic)
are generated.

..  index:: pair:Device;simulation

Step 4: Run TiberCAD
--------------------

Now we can run TiberCAD:

  ``tibercad bulk.tib``

The generated Output files are:

:file:`driftdiffusion_materials.dat` : material (mesh) regions, in this case just region 1

:file:`driftdiffusion_nodal.dat` : nodal quantities (here conduction and valence band)

:file:`sweep_driftdiffusion_Vb.dat` : integrated current at the two contacts for each sweep step.


.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _Tutorial_00: http://www.tibercad.org/documentation/tutorial/tutorial_00_bulk_si

.. rubric:: Footnotes

