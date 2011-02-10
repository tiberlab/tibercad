.. _GettingStarted2D:

Getting started 2D
==================


In this second example we will refer to the :ref:`Tutorial 04` (Si n-Mosfet) that you can find
in the :ref:`Tutorial` directory.

.. _Step2D: 

Step 1: Modeling the device
---------------------------

Again, as a first step, we have to model the device.
   We'll see in some details how to design and mesh a mosfet device with GMSH.

      * In the GMSH script ``mosfet.geo`` , several variables are defined and given a value in this way:
   ::

     lsub=0.03;
     lacc=0.002;
     lct=0.0005;
     lg=0.0015;
     lh=0.01;
     lc=0.0005;

these variables are used in the script to assign proper values to the **mesh characteristic length** of the defined Points::

  Lg_2 = 0.0375;
  d = 0.01;
  Ls = 0.1;
  h = 0.25;
  b = 0.0025;
  o = 0.005;
  xd = Lg_2 + d;
  xd2 = Lg_2 + d / 2;
  xmax = xd + Ls - d;
  
These other convenient variables are used to parametrize the most relevant geometrical
features, such as channel length, oxide thickness, and so on.

  * Geometrical **Points** and **Lines** are defined to design the device structure; the
  fourth parameter in **Point assignment** is the characteristic length associated to
  that point: this is an essential feature to control the mesh density and refine it
  where necessary (usually n the channel region).

.. note:: |warn|

   In a 2D simulation it is assumed that the geometrical model is restricted to the
   **xy-plane (z = 0)**. Any other geometrical orientation could give unpredictable results.
    
::
   
  Point(1) = {0, -h, 0, lsub};
  Point(2) = {0, 0, 0, lc};
  Point(3) = {xmax,-h,0.0,lsub};
  Point(4) = {-xmax,-h,0.0,lsub};
  Point(5) = {xmax,0,0.0,lh};
  Point(6) = {-xmax,0,0.0,lh};
  ..........................
  Line(1) = {4,1};
  Line(2) = {3,13};
  Line(6) = {4,14};
  Line(7) = {10,9};
  Line(8) = {12,2};
  Line(9) = {8,7};
  Line(10) = {11,8};
  Line(11) = {9,12};
  Line(13) = {7,6};
  ..........................
  
----

  * Definition of a surface: first a **line loop** is composed, listing all the lines constituting 
  the boundary of the surface; then this line loop is assigned to a **Plane
  Surface** object (this procedure can be alternatively performed through the graphical
  interface)::
  
    Line Loop(40) = {28,2,-34,33,8,29,-31,-30,-6,1};
    Plane Surface(41) = {40};
    ..........................
    
  * Definition of the **Physical surfaces** : each of them is composed by one or more
  geometrical **Plane Surface**. For example, **Physical surface 2** comprises the
  two separated contact regions, while **Physical surface 3** corresponds to the oxide
  region.

The Physical surfaces are the 2D Physical regions of the mesh and will be as
signed to the related TIBERCAD regions through the keyword *mesh_regions* ( see :ref:`Step2D` )::

    Physical Surface(1) = {41}; // n-Si
    Physical Surface(2) = {44,47}; // n+-Si
    Physical Surface(3) = {46}; // SiO2
  
* Definition of the **Physical Lines**: In this 2D simulation, 1D physical regions are
used to carry information about boundary condition regions. In other word, each
**Physical Line** corresponds to a boundary condition (a contact in the case of a
driftdiffusion calculation): thus Physical Line 1 refers to source contact, P.L. 2 to
gate contact, P.L. 3 to drain contact. The numerical identifications of these **Physical Lines** 
will be assigned to TIBERCAD **BC regions** by means of the *BC_reg_numb* instruction::
  
    Physical Line(1) = {13}; // source
    Physical Line(2) = {39,38}; // gate
    Physical Line(3) = {19}; // drain
    
In :ref:`fig.3.1` the obtained geometrical model is shown.

Step 2: Meshing the device
--------------------------

The *.geo* script file with the geometrical description can be run in GMSH, to display
the modelled device and to mesh it through the GMSH graphical interface ( see :ref:`fig.3.2` ).
Alternatively, a *non-interactive* mode is also available in GMSH, without graphical user
interface. 

For example, to mesh this 2D tutorial in non-interactive mode, just type::

  gmsh mosfet.geo -2 -o mosfet.msh

  
Step 3: TiberCAD Input file
---------------------------

Now we have to write down the TIBERCAD **input file** (see *mosfet.tib* in the Tutorials).

.. _fig.3.1: 

.. figure:: geomosfet.*
   :align: center
   
   fig. 3.1 (Geometrical structure as defined by GMSH modeller)

.. _fig.3.2: 

.. figure::meshmosfet.*
   :align: center

   fig. 3.2 (2D Mesh for the mosfet device obtained with GMSH)


1 - Definition of Device Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three TIBERCAD regions are defined: to each of them, one mesh region is associated,
that is the Physical Surfaces **1**, **2** and **3** defined in :ref:`Step2D` . However, in general more than
one mesh region can be associate to a single TIBERCAD region, if this is convenient::

  Region substrate
    {
     mesh_regions = 1
     material = Si
     doping = 1e18 doping_type = acceptor
    }
  Region contact
    {
     mesh_regions = 2
     material = Si
     doping = 5e19 doping_type = donor
    }
  Region oxide
    {
     mesh_regions = 3
     material = SiO2
    }


2 - Definition of Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we define the **Simulation** *dd*: it belongs to the class **driftdiffusion**::

  model driftdiffusion
    {
     options
       {
        simulation_name = dd
        physical_regions = all
       }
  
We declare two **driftdiffusion physical models**: the first defines a *srh* recombination model ( see :ref:`srh` ); 
the second defines a *field-dependent* **mobility** model for
electrons which implements a *doping dependence* for the **low-field mobility** ( see :ref:`mobility` )::

  physical_model recombination
    {
     model = srh
    }
  physical_model electron_mobility
    {
     model = field_dependent
     low_field_model = doping_dependent
    }
    
  
3 - Definition of Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source, drain and gate contacts of the Mosfet device are defined as **Boundary
conditions regions** (*BC_ Region source* , *BC_Region drain*, *BC_Region gate*) in the
following way::

  BC_Region gate
    {
     BC_reg_numb = 2
     type = schottky
     barrier_height = 3.0
     voltage = @Vg[0.0]
    }
  BC_Region source
    {
     BC_reg_numb = 1
     type = ohmic23
     voltage = 0.0
    }
  BC_Region drain
    {
    BC_reg_numb = 3
    type = ohmic
    voltage = @Vd[0.5]
    }

To each of the **BC regions**, one *BC_reg_numb* is assigned, that is one of the Physical
Lines **1,2,3** defined in  :ref:`Step2D` , which represent the contact regions.

.. note:: |warn|

   Note that, while *source* and *drain* are defined as *type = ohmic, gate* **BC_region** is
   defined as *type = schottky; barrier height = 3.0* specifies the metal/oxide interface barrier
   and depends on the contact metal workfunction.

Drain voltage is defined as *@Vd[0.5]* and gate voltage as *@Vg[0.0]*. This specifies
that the value of the voltage will be determined at each moment of the simulation, by
the value of the two variables :math:`V_d` and :math:`V_g`, which will be assigned in the **sweep** definition.

4 - Definition of Simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two sweeps are requested for this simulation, that is an external loop on :math:`V_g` (the gate
voltage) and an internal loop on :math:`V_d` (the drain voltage) for each value of :math:`V_g`; in this way,
the IV drain characteristics for a series of gate biases are obtained in output::

  sweep_1
    {
     simulation = dd
     variable = Vd
     start = 0.0
     stop = 2.0 #0.1
     steps = 200 #200 #1
     # plot_data = true
    }
  
  sweep_2
    {
     variable = Vg
     start = -0.1
     stop = 0.5
     steps = 6
     simulation = sweep_1
     #simulation = driftdiffusion
    }
  
5 - Definition of Execution parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the **Simulation** section , we decide the simulation dimension (*dimension = 2*), then
which simulations to perform and in which order ; we set *solve = sweep_2*, to execute the
external gate voltage sweep *sweep_2* which in its turn call the sweep *sweep_1* where drain
current is calculated for all the chosen drain voltage steps by running *dd* **simulation**::

  $Simulation
    {
     meshfile = mosfet.msh
     dimension = 2
     temperature = 300
     solve = sweep_2
     resultpath = output_IV_char
     output_format = vtk
     plot = (Ec, Ev, QFermi_e, QFermi_h, eDensity, hDensity, eCurrent, hCurrent,
     NetRecombination, EField, ElPotential, ContactCurrents)
    }
  
Output files with conduction and valence band profiles, quasi-fermi levels, electron
and hole density, recombination, electric field and potential (*plot = :math:`E_c,E_v` ,.....*) will be
generated, together (*ContactCurrents*) with a file with all the calculated values of the
drain current at the contacts for each gate bias step (the IV characteristics).


Step 4: Run TiberCAD
--------------------

Now we can run TiberCAD:

  > tibercad mosfet.tib

The generated Output files are:
  
:file:`driftdiffusion_materials.vtk`: information about the material regions of the device.

:file:`driftdiffusion_nodal.vtk` : output for the nodal quantities which have been calculated, 
e.g. conduction and valence bands, (quasi)fermi levels, electron and hole density and mobility.

:file:`driftdiffusion_elemental.vtk` : output for the elemental quantities (e.g. electric field, 
current density).

:file:`sweep_2_driftdiffusion_Vg_0.000_Vd.dat` and similar for all the :math:`V_g` steps: drain
current characteristics for each :math:`V_g` bias.


.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%
          
.. _Tutorial 04: http://www.tibercad.org/documentation/tutorial/tutorial_04_n_mosfet_2d

.. rubric:: Footnotes