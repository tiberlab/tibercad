..   <marker>

.. _DriftDiffusionGetting:

Drift Diffusion
=================================================

Introduction
-------------

The module Drift-Diffusion simulates electron and hole transport based on the drift-diffusion
approximation (see [Selberherr]_ ).

The system of equations that is solved is given by

.. math::
   :label: dd_eq_ddsystem
   
   -\nabla(\varepsilon\nabla\varphi - \mathbf{P}) & =  -e(n - p - N_d^+ + N_a^-) \\
   -\nabla(\mu_n n ( \nabla\phi_n + P_n \nabla T)  ) & =  R \\
   -\nabla(\mu_p p (\nabla\phi_p + P_p \nabla T) ) & =  -R, 

where :math:`P` is the electric polarization due to e.g. piezoelectric effects and :math:`R` is the net 
recombination rate, i.e. recombination rate minus generation rate. :math:`P_n` and :math:`P_p` are the electron
and hole thermoelectric powers, respectively. The models for the mobilities and the net
recombination rates can be specified in the ``Physics`` section.
Recombination and trap models can be specified both for bulk and interface device regions.



Example
--------

This very simple example presents the basic setup of a Drift-Diffusion simulation
by simulating the current in a piece of bulk silicon.

The following files should be in your working directory:

  bulk.tib_  :            tiberCAD input file 

  bulk.msh_ :             mesh file

The mesh file can be obtained from the following ``GMSH`` script (bulk.geo_):

::

   L = 1;
   d = 0.01;

   Point(1) = {0, 0, 0, d};
   Point(2) = {L, 0, 0, d};

   Line(1) = {1, 2};

   Physical Line("bulk") = {1};
   Physical Point("anode") = {1};
   Physical Point("cathode") = {2};


The simulated device is a homogeneous piece of bulk silicon with a length of :math:`1\, \mu m`.
At the two ends we define the contacts (``anode`` and ``cathode``).

The tiberCAD input file (bulk.tib_) is shown in the following listing:

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

  Module driftdiffusion
  { 
    plot = (Ec, Ev, eQFermi, hQFermi, ContactCurrent)

    Contact anode { voltage = $Vb }
    Contact cathode { }
  }

  Module sweep
  {
    solve = driftdiffusion
    variable = $Vb
    start = 0.0
    stop = 1
    steps = 10
    plot_data = true
  }

  Simulation
  {
    verbose = 2

    solve = sweep

    resultpath = output
    output_format = grace
  }


In the ``Device`` block we provide the name of the mesh file (``bulk.msh``), assuming default mesh units
of micrometers (``mesh_units = 1e-6``).
The model has a single region ``bulk`` with material silicon (``material = Si``), which is slightly n-doped.

Next we specify the simulation model. Here we only solve the Drift-Diffusion model (``Module = driftdiffusion``).
We do not specify any physical models for recombination or mobility in this example, but we have to define the electrical contacts
using the ``Contact`` keyword.
At the anode we apply a variable voltage by providing a variable name (``$Vb``).
We also specify which quantities to plot.

In this simulation we want to calculate the IV characteristic of a piece of bulk Si.
For this porpose we define a sweep on the anode voltage (``Module sweep``).
In the sweep block we have to specify the sweep variable, which in our case is ``$Vb``.
The ``plot_data = true`` option will make the sweep plot all results after each sweep step.

Finally, in the ``Simulation`` block we define the output directory and format (in this case ``grace`` for 1D ASCII files),
and we specify that we want to solve the sweep.


The generated Output files are:


* ``driftdiffusion_msh.dat``      : mesh-dependent quantities (conduction and valence band edges and quasi Fermi levels)

* ``driftdiffusion_Vb_<step>_msh.dat``: mesh-dependent quantities for each sweep step 
 
* ``sweep_driftdiffusion.dat``   : integrated current at the two contacts for each sweep step



..  _bulk.geo: http://www.tiberlab.com/images/stories/products/device_sim/tibercad/manuals/resources/bulk_0.geo
..  _bulk.msh: http://www.tiberlab.com/images/stories/products/device_sim/tibercad/manuals/resources/bulk_0.msh
..  _bulk.tib: http://www.tiberlab.com/images/stories/products/device_sim/tibercad/manuals/resources/bulk_2.tib




 


.. rubric:: Footnotes


..   </marker>


