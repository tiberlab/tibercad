..   <marker>

.. _DscGetting:


Simulation Dye Solar Cells
==============================

Model DSC
^^^^^^^^^

The DSC model is tagged as ``dssc`` . In **options** subsection we indicate the simulation name::

  model dssc
    {
     options
       {
        simulation_name = <name_of_the_model>
        plot(Potential, Density, Current, ContactCurrents)
       }
     BC_regions
       {
        ...
       }
    }
    
Boundary Conditions
^^^^^^^^^^^^^^^^^^^

The boundary conditions are inserted in the **Model** section, the two contacts are the
photoanode and the cathode. The photoanode must be the boundary of a region where
:math:`\bf{TiO_2}` is present, on the contrary the cathode must be the boundary of a region where
the **electrolyte** is present.

::

  BC_regions
    {
     BC_region <name_of_the_anode>
       {
        type = ohmic
        bias = @V[0.0]
       }
     BC_region <name_of_the_cathode>
       {
        type = Pt
        load = @R[1e8]
       }
    }
    
The anode is an ohmic contact, it contains only the sweep of the bias applied for the
electrochemical potential of the electrons. The cathode must contain the initial external
resistance ( **load = R** ). If a simulation of the cell under illumination is made the initial
value of *load* must be set to a high value, in the range of :math:`10^8 \omega cm^2` in order to have 
no current during the light intensity sweep. In case instead the simulation performed is
under dark conditions where an external bias only is applied ``load = 1.0`` .

Generation model
^^^^^^^^^^^^^^^^

The generation term is related to the flux of photons which reaches the active :math:`TiO_2` 
regions and the dye present in the cell. We assume a simple Beer-Lambert exponential
decay for charge generation of the form:

..  math::
    :nowrap:
    :label:
    
    \begin{equation}\label{generation}
    G(\mathbf{r}) = \int^{\lambda_{2}}_{\lambda_{1}} \alpha(\lambda)
    \Phi (\lambda) e^{-\alpha(\lambda) \hat{n} \cdot \mathbf{r}}
    d\lambda
    \end{equation}

where :math:`\alpha` is the absorption coefficient (in :math:`\mu^{-1}` ) 
of the chosen Dye, :math:`\Phi(\lambda)` the intensity of the light power at that 
wavelength for the spectrum of the sun.
The part of the input file for the generation model::

  model dssc_generation
    {
     options
       {
        regions = (<TiO2_region_name_1>, <TiO2_region_name_2>, ...)
        plot = (Distance, Generation)
        light_direction = <vector_indicating_the_direction_of_light>
       (example, illumination from x direction: light_direction = (1, 0, 0))
        light_intensity = @x
        dye = N719
       }
     BC_Regions
       {
        BC_Region <name_of_the_boundary>
          {
          }
       }
    }
    
In the generation input flle must be specified:

*  ``regions`` : the regions where we want that there is generation (where the Dye is present);

*  ``plot`` : if we want to plot the generation;

*  ``light_direction`` : the vector which fixes the direction from where the light comes;

*  ``light_intensity`` : the light intensity;

*  ``light_intensity`` : the light intensity;

*  ``dye`` : the dye used in the cell;


The light intensity is defined in a sweep ( see section :ref:`Solver_Dsc` ). It can be set to reach 1
that means one Sun, or a larger or smaller illumination intensity (0.1, 2.0, etc.).
There is another flag called ``illumination_spectrum`` that can be used if we want to
change the spectrum profile F with another illumination source (for example if we assume
the cell under concentration of light). The file used by default for :math:`\Phi` is the standard 1.5
AM spectrum of the sun contained in the material database in the file ``Sun1p5am`` .

The last part of the generation is the definition of the boundary. This boundary says
to the code which is the boundary region of the grid from where the light penetrates in
the device. The information given by the boundary and the light direction vector defines
the coming of light and direction of rays.

Device
------

Device section for a DSC simulation::

  Device
    {
     Region <name of the porous region>
       {
        mat = TiO2mes
        porosity = 0.5
       }
     Region <name of the electrolyte region>
       {
        mat = TiO2mes
        TiO2 = false
       }
    }
    
For every region of the device it is specified the material file from the database (the
TiO2mes contains standard parameters for both :math:`TiO_2` and electrolyte, so it can be used
for both regions). For every region must be specified if it contains :math:`TiO_2`, or electrolyte
or both. 

This can be done setting two flags called ``TiO2`` and ``electrolyte`` . If set true
the material is present in the region, if set **false** it is not present. By default they are
assumed both **true** (porous region). In the second region of the example showed there
is electrolyte only and ``TiO2`` = **false** is explicitly specified. In case both materials are
present (porous region) a porosity must be specified (in the range between 0, :math:`TiO_2` only,
and 1, electrolyte only). 

If one material is not present the porosity is automatically set
to 0 or 1.

Physics
-----------

The set of parameters that can be defined for the entire device are enlisted in table 3.1.
For the generation::

  dssc
    {
     generation = dssc_generation
    }
    
.. _Solver_Dsc:

Solver
----------

Three sweeps are needed for the plot of the entire I-V under illumination. The first sweep
is needed to make the transition from dark conditions to full open-circuit conditions under

.. _dsc_parameters:

..  math::
    :nowrap:
    :label: 
    
    \begin{table}[!htb]
    \center
    \begin{tabular}{l|p{8cm}|l}
    \multicolumn{3}{c}{\textbf{Parameters}} \\
    \hline
     & \textbf{Poisson equation} & \\
    \hline
    \texttt{porosity} & Porosity of porous material &  \\
    \texttt{perm\_oxide} & Relative perm. of the TiO$_2$ &  \\
    \texttt{perm\_electrolyte} & Relative perm. of the electrolyte & \\
    \texttt{k\_3} & Oxidized Dye regeneration rate constant & s$^{-1}$ \\
    \hline
     & \textbf{Drift Diffusion equation} & \\
    \hline
    \texttt{ne} & Dark electron density & cm$^{-3}$ \\
    \texttt{nI} & Dark iodide density & mol/l \\
    \texttt{nI3} & Dark electron density & mol/l \\
    \texttt{mu\_e} & Electron mobility & cm$^{2}$/Vs \\
    \texttt{D\_I} & Iodide diffusion constant & cm$^{2}$/s \\
    \texttt{D\_I3} & Triiodide diffusion constant & cm$^{2}$/s \\
    \texttt{D\_C} & Cation diffusion constant & cm$^{2}$/s \\
    \hline
     & \textbf{Recombination term} & \\
    \hline
    \texttt{k\_e} & Recombination constant rate & s$^{-1}$ \\
    \texttt{beta} & Non-linearity exponent in the electron density recombination &  \\
    \end{tabular}
    \caption{Physical parameters that can be set for the model divided
    in subsets relative to different processes in the cell.}
    \label{table:dsc_parameters}
    \end{table}


illumination. The high value of the load R maintains the current to zero. Then a second
sweep is used to pass from open circuit condition to short circuit condition lowering the
external load from a high external load to a small one (R = 1). Finally, the voltage
sweep compute the I-V characteristics under illumination. In case of dark simulation
(application of an external bias without illumination) the first and second sweep are not
needed. The external load for the cathode can be set directly to R = 1.

::

  Solver 
    {
     Sweep
       {
        sweep_gen
          {
           simulation = (dssc_generation, dssc)
           variable = x
           values = (0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1)
           plot_data = true
          }
        sweep_R
          {
           simulation = dssc
           variable = R
           values = ( 1000, 100, 1)
           plot_data = true
           }
        sweep_V
          {
           simulation = dssc
           variable = V
           values = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0,8, 0.9, 1.0)
           plot_data = true
          }
        }
    }
    
The intensity of illumination can be changed sweeping the value of x different to 1 (1 =
1 Sun of power).

Output
----------

The output that want to be plotted are enlisted in section Model within option. See 
:ref:`Table nodal<dsc_nodal>` , :ref:`Table elemental<dsc_elemental>` and :ref:`Table scalar<dsc_scalar>` .

..  _dsc_nodal:

..  math::
    :nowrap:
    :label:
    
    \begin{table}[!htb]
    \center
    \begin{tabular}{l|p{8cm}|l}
    \multicolumn{3}{c}{\textbf{Nodal quantities}} \\
    \hline
     & \texttt{Potential} & \\
    \hline
    \texttt{eQFermi} & Electron electrochemical potential & eV ($-e\phi_n$)\\
    \texttt{IQFermi} & Iodide electrochemical potential & eV ($-e\phi_{I^{-}}$)\\
    \texttt{I3QFermi} & Triiodide electrochemical potential & eV ($-e\phi_{I^{-}_{3}}$)\\
    \texttt{CQFermi} & Cation electrochemical potential & eV ($-e\phi_C$)\\
    \texttt{ElPotential} & Electrostatic potential & eV \\
    \texttt{Eredox} & Electrolyte electrochemical potential & eV ($-e\phi_{Red}$) \\
    \hline
     & \texttt{Density} & \\
    \hline
    \texttt{eDensity} & Electron density & cm$^{-3}$ \\
    \texttt{IDensity} & Iodide density & cm$^{-3}$ \\
    \texttt{I3Density} & Triiodide density & cm$^{-3}$ \\
    \texttt{CDensity} & Cation density & cm$^{-3}$ \\
    \texttt{Generation} & The net electron generation rate & cm$^{-3}$s$^{-1}$ \\
    \texttt{NetRecombination} & The net recombination rate & cm$^{-3}$s$^{-1}$ \\
    \end{tabular}
    \caption{Nodal quantities. The flags \texttt{Potential} and
    \texttt{Density} allow to plot all the electrochemical potentials
    and all the densities, including generation and recombination.}
    \label{table:dsc_nodal}
    \end{table}


.. _dsc_elemental:

..  math::
    :nowrap:
    :label: 

    \begin{table}[!htb]
    \center
    \begin{tabular}{l|p{5cm}|l}
    \multicolumn{3}{c}{\textbf{Elemental quantities}} \\
    \hline
     & \texttt{Current} & \\
    \hline
    \texttt{ElField} & Electric Field & Vcm$^{-1}$ \\
    \texttt{CurrentDensity} & Total current density & Acm$^{-2}$ \\
    \texttt{eCurrentDensity} & Electron current density & Acm$^{-2}$ \\
    \texttt{ICurrentDensity} & Iodide current density & Acm$^{-2}$ \\
    \texttt{I3CurrentDensity} & Triiodide current density & Acm$^{-2}$ \\
    \texttt{CCurrentDensity} & Cation current density & Acm$^{-2}$ \\
    \end{tabular}
    \caption{Elemental quantities. The flag \texttt{Current} allows to
    plot all of them in the x, y and z components.}
    \label{table:dsc_elemental}
    \end{table}

..  _dsc_scalar:

..  math::
    :nowrap:
    :label: 

    \begin{table}[!htb]
    \center
    \begin{minipage}{8cm}
    \center
    \begin{tabular}{l|l|l}
    \multicolumn{3}{c}{\textbf{Scalar quantities}} \\
    \hline
    \texttt{ContactCurrents} & Contact currents & *\footnote{depends on dimension and symmetry}\\
    \end{tabular}
    \end{minipage}
    \caption{Scalar quantities.} \label{table:dsc_scalar}
    \end{table}


.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%

..   </marker>
