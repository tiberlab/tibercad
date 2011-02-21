.. _DriftDiffusion:

Drift Diffusion
=================================================


In this example we will see a very simple TiberCAD simulation:

1D calculation of Poisson and drift-diffusion for a bulk Silicon sample.

The following files should be in your working directory:

|   bulk.tib_ : **TiberCAD** input file 
|
|   bulk.msh_ : mesh file

The mesh file can be obtained from the following GMSH geo file : bulk.geo_ 

This is the summary of the Sections of this Tutorial :

*  :ref:`tut0step1` 
*  :ref:`tut0step2` 
*  :ref:`tut0step3` 
*  :ref:`tut0step4` 
*  :ref:`tut0step5` 
 

..  _tut0step1:
 
Step 1 - Modeling the device
----------------------------

As a first step, we have to model the device. To do so, you can use DEVISE 
module of ISE-TCAD 9.5 software package or GMSH program.

Here we'll see in details the procedure for GMSH.

There are two possible ways to use GMSH:

Interactive, using the graphical interface
Using a script file
In the following we'll see how to write a basic GMSH script ( bulk.geo_ ); 
for any details please refer to GMSH manual GMSH (http://geuz.org/gmsh/).

 
| 
| 


  |warn|: In a **1D** simulation it is assumed that the geometrical model is restricted to the **x axis** .

          In a **2D** simulation it is assumed that the geometrical model is restricted to the **xy-plane (z=0)** 

          Any other geometrical orientation could give impredictable results

 

In a GMSH script, several variables can be defined and given a value in this way::

  L = 1
  d = 0.01
  
these are valid GMSH variables: L is just the length of the Si sample; d is the value 
of a *characteristic mesh length* (see below).

Definition of geometrical entities **Points** ::

  Point(1) = {0, 0, 0, d}
  Point(2) = {L, 0, 0, d}

The first three expressions inside the braces on the right hand side give the three 
X, Y and Z **coordinates** of the point; the last expression **(d)** sets the **characteristic 
mesh length** at that point, that is the *"size"* of a mesh element, defined as the length 
of the segment for a line segment, the radius of the circumscribed circle for a triangle 
and the radius of the circumscribed sphere for a tetrahedron.

Thus, the smaller is the value of d, the greater is the mesh density close to that point. 

The size of the mesh elements will then be computed in GMSH by linearly interpolating 
these characteristic lengths in the whole mesh.
 

Definition of a geometrical entity **Line** ::

  Line(1) = {1, 2}

The two expressions inside the braces on the right hand side give the identification 
numbers of the start and end points of the line.


..  figure:: ../data/bulk_geo_fig.png
    :align: center
    :scale: 50%


Definition of the physical entity **Physical Line bulk** ::

  Physical Line("bulk") = {1}

The expression(s) inside the braces on the right hand side give the identification numbers 
of all the **geometrical lines** that need to be grouped inside the *physical line* .

In this way, in general, physical regions are created which associate together geometrical regions, 
and then the related mesh elements, which share some common physical properties. It's only these 
physical regions which can be referred to outside GMSH. In TiberCAD, this is done by associating 
one or more physical regions to a **TiberCAD region** through the keyword *mesh_regions* (see in the following).

Definition of two physical entities Physical Point::

  Physical Point("anode2)   = {1}
  Physical Point("cathode") = {2}
 

  |note| : In general, in a nD simulation, **(n-1)D** physical regions (points in 1D, lines in 2D, surfaces in 3D) 
   are used by TiberCAD to impose the required boundary conditions.

   Each (n-1)D physical region defined in this way in GMSH will be associated in TiberCAD to a boundary condition region, 
   through the keyword **BC_reg_numb** . Thus, in this case, Physical points Anode and Cathode will be associated respectively 
   to two **Contact** regions (see in the following).
 

..  _tut0step2:

Step 2 - Meshing the device 
---------------------------
 

The *.geo* script file with the geometrical description can be run in GMSH, to display the modelled device 
and to mesh it through the GMSH graphical interface.
Alternatively, a *non-interactive* mode is also available in GMSH, without graphical user interface. 

For example, to mesh this 1D tutorial in non-interactive mode, just type:

  |idea|    gmsh bulk.geo -1 -o bulk.msh

where bulk.geo_ is the geometrical description of the device with GMSH syntax:

  -1 means *1D mesh generation*

some command line options are:

  -1 , -2, -3 to perform *1D, 2D or 3D* mesh generation,

  -o **mesh_file.msh** to specify the name of the mesh file to be generated

In this way, a .msh has been generated and is ready to be read in TiberCAD.

..  _tut0step3:

Step 3 - TiberCAD Input file
-----------------------------


Now we have to write down the **TiberCAD input file** ( bulk.tib_ ). For a detailed reference see the user guide 
(Input file) and Getting started 1

Description of Device Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we have to list all the **TiberCAD Regions** present in our Device: a TiberCAD Region is usually a section 
of the device featuring the same material and possibly the same doping.

::

  Device
    {
     meshfile = bulk.msh

     Region bulk 
       {
        material = Si

        Doping
          {
           Nd = 1e16
           type = donor
          }
       }
    }
    
The TiberCAD Region bulk is made of Silicon and n-doped with a concentration :math:`1 x 10^{16} cm^{-3}` .

Through the keyword **Region** , one GMSH physical region (Physical Lines in 1D, Physical Surfaces in 2D, 
Physical Volumes in 3D) previously defined in the GMSH mesh ( :ref:`tut0step1` ), can be associated 
to the present TiberCAD Region, in this way::


  Region  GMSH_physical_region_name


In this case, the Physical Line bulk is associated to the TiberCAD Region bulk.

Alternatively, through the optional keyword **mesh_regions** , one or more GMSH physical regions can be 
associated to a single TiberCAD Region .

 

Definition of Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we define the **Simulation** driftdiffusion_1: it belongs to the class **driftdiffusion** 

::

  Module driftdiffusion
    {               
     name = driftdiffusion  # this is the default name

     regions = all          # 'all' is the default

     # what we want to plot
     plot = (Ec, Ev, eQFermi, hQFermi, ContactCurrent)
     ....
       
The TiberCAD simulation *driftdiffusion_1*  , belonging to the model driftdiffusion, will be applied 
to the whole device structure (``regions = all``)

 

Definition of Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The anode and cathode contacts of our 1D Si sample are defined as **Boundary conditions regions** 
(``Contact anode, Contact cathode``) in the following way::


Module driftdiffusion
{ 
              
  name = driftdiffusion  

  regions = all

  plot = (Ec, Ev, eQFermi, hQFermi, ContactCurrent)

  Contact anode { voltage = $Vb }
  Contact cathode { }

}
    
Through the keyword Contact , one (n-1) -dimension GMSH physical region (Physical Point in 1D, 
Physical Line in 2D, Physical Surface in 3D) previously defined in the GMSH mesh ( :ref:`tut0step1` ), can be 
associated to the present TiberCAD Contact, in this way:

Contact  GMSH_physical_region_name 

In this case, the *Physical Point* **anode** is associated to the TiberCAD Contact anode and the Physical 
Point cathode is associated to the TiberCAD Contact cathode.

Both contacts are defined as *ohmic* , cathod is assigned a fixed ``voltage = 0.0`` , while anode voltage is given 
by the value of the variable *Vb* ::


  voltage = @Vb


Definition of Simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vb is specified in the sweep block, in the Solver section ::


Module sweep
{
  solve = driftdiffusion
  variable = $Vb
  start = 0.0
  stop = 1
  steps = 10
  plot_data = true
}
    
In this way, the simulation driftdiffusion_1 is performed for 10 ( steps = 10) values of the anode 
voltage (variable = Vb), between 0 and 1.

For each step we want to plot the solution variables specified in the driftdiffusion module 
(plot_data = true).
 

Definition of Execution parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In the **Simulation** section , we decide *which* simulations to perform and in which *order*; set ``solve = sweep`` , 
to execute the sweep which run driftdiffusion_1 simulation for the specified loop.

::

  Simulation
    {
     verbose = 2 

     solve = sweep

     resultpath = output
     output_format = grace
    }

Output files with conduction and valence band profiles (plot = Ec,Ev..) and all the calculated values 
of the current at the contacts (*ContactCurrents*) (the IV characteristic) are generated.

To increases the amount of information written to the screen we can vary the verbose level (verbose = 2).
 
..  _tut0step4:

Step 4 - Run TiberCAD
-----------------------

Now we can run TiberCAD

|  by double clicking on bulk.tib file (in Windows)
|  
|  or by command line in linux: tibercad bulk.tib


..  _tut0step5:

Output 
----------

The generated Output files are:

* driftdiffusion_materials.dat  : material (mesh) regions, in this case just region 1

* driftdiffusion_nodal.dat      : nodal quantities (here conduction and valence band)

* sweep_driftdiffusion_Vb.dat   : integrated current at the two contacts for each sweep step
 


Attachment	Size

* bulk.tib_	1.16 KB
* bulk.geo_	181 bytes
* bulk.msh_	4.49 KB

 


 
..  _bulk.tib: http://www.tibercad.org/files/bulk_2.tib
..  _bulk.geo: http://www.tibercad.org/files/bulk_0.geo
..  _bulk.msh: http://www.tibercad.org/files/bulk_0.msh

.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%

.. rubric:: Footnotes
