
.. _DscTheory:


Simulation Dye Solar Cells
========================================


Introduction
------------------------------

For a brief list of literature to understand Dye Solar cells (DSC) see [Kalyanasundaram]_ . 

For a review of the model see [Gagliardi]_ .

The model consists in a set of drift-diffusion equations for the propagation of ions
and electrons coupled with Poisson equation:

.. math::
   
   :label:

   \begin{eqnarray}
    \nabla \cdot (\mu_{e}n_{e}\nabla\phi_{e}) & = & (G - R) \\
    \nabla \cdot (\mu_{I^{-}}n_{I^{-}}\nabla\phi_{I^{-}}) & = & -\frac{3}{2}(G - R) \\
    \nabla \cdot (\mu_{I^{-}_{3}}n_{I^{-}_{3}}\nabla\phi_{I^{-}_{3}}) & = & \frac{1}{2}(G - R) \\
    \nabla \cdot (\mu_{c}n_{c}\nabla\phi_{c}) & = & 0,
   \label{eq1}
   \end{eqnarray}

where :math:`\mu_{\alpha}` refers to carrier mobilities, :math:`n_{\alpha}` to charge concentrations and :math:`\phi_{\alpha}` to 
electrochemical potentials. **R** is the recombination term and **G** the generation term due to illumination. 
In order to take into account the trap density we use a density dependent
mobility developed from multi-trapping model:

.. math::
   
   :label:

   \begin{equation}\label{difftraps}
    \mu_{e}(n_{e}) = \mu_{0} \left ( \frac{n_{e}}{N_{t}} \right )^{\frac{1 - a}{a}}
   \end{equation}

|

where :math:`a` is the trap exponent, :math:`N_t` the trap density and :math:`\mu_0` a constant. The energy trap
density is assumed to form an exponential tail below the conduction band of the semiconductor:

.. math::
   
   :label:
   
   \begin{equation}\label{denstrap}
    g_{T}(E) = \frac{a N_t}{kT} e^{\frac{-a E}{kT}}.
   \end{equation}

|

The Poisson equation to handle the internal potential drop:

.. math::
   
   :label:
   
   \begin{equation}\label{poisson}
    -\varepsilon\triangle \varphi =  q[n_{c} + N_{D}^{+} - n_{I^{-}} - n_{I_{3}^{-}} - (n_{e} - \bar{n}_{e})],
   \end{equation}

|

where N+ D is the amount of ionized dyes and it is equal to:

.. math::
   
   :label:
   
   \begin{equation}\label{dyeion}
    N_{D}^{+} = \frac{G}{k_{3}}
   \end{equation}

|

with **G** the generation term and :math:`k_{3}` the rate constant of dye regeneration. The dielectric
constant, :math:`\varepsilon` , of the mesoporous material is a mix the two dielectric functions of the semiconductor and the electrolyte. 
We use the Maxwell-Garnet model where the dielectric
function of the mixed medium becomes:

.. math::
   
   :label:
   
   \begin{equation}\label{diel}
    \varepsilon = \varepsilon_{s}\frac{\varepsilon_{e} + 2\varepsilon_{s} + 2\epsilon_{p}\varepsilon_{e} - 2\varepsilon_{s}\epsilon_{p}}
    {\varepsilon_{e} + 2\varepsilon_{s} -\epsilon_{p}\varepsilon_{e} +\varepsilon_{s}\epsilon_{p}}
   \end{equation}

|

where :math:`\varepsilon_{s}` and :math:`\varepsilon_{e}` are the dielectric constants of the semiconductor and the electrolyte,
respectively, and :math:`\epsilon_{p}` is the porosity of the medium. The recombination term depends
largely on the loss mechanisms at the electrolyte/oxide interface which follows the reaction chain:

.. math::
   
   :label:
   
   \begin{eqnarray}
   I^{-} & \rightleftharpoons & I + e \\
   2I & \rightleftharpoons & I_{2} \\
   I_{2} + I^{-} & \rightleftharpoons & I^{-}_{3}.
   \label{reaction_loss}
   \end{eqnarray}

|

From the chemical path we can get a formula for the interface recombination (considering
that the first chemical reaction is the slow process):

.. math::
   
   :label:
   
   \begin{equation}\label{ricombinazione}
    R = k_{e} \left [  \left ( \frac{n_{e}}{\bar{n}_{e}} \right )^{\beta}\bar{n}_{e}\sqrt{\frac{n_{I^{-}_{3}}}{n_{I^{-}}}}
    - \bar{n}_{e}\sqrt{\frac{\bar{n}_{I^{-}_{3}}}{\bar{n}^{3}_{I^{-}}}} n_{I^{-}}\right
    ],
   \end{equation}

|

where the electron rate :math:`k_{e}` is the recombination rate constant.


For the boundary conditions of the model we assume at the photoanode:

*  :math:`\phi_{e} = V`: electrochemical potential of electrons set with the voltage applied;
*  :math:`\nabla\phi_{I^{-}} = 0`: no iodide current at the photoanode;
*  :math:`\nabla\phi_{I_{3}^{-}} = 0`: no triiodide current at the photoande;
*  :math:`\nabla\phi_{c} = 0`: no cationic current;
*  :math:`\nabla\varphi = 0`: no charged layer at the photoanode;


at the cathode:

*  :math:`\nabla\phi_{e} = 0`: no electronic current at the cathode;
*  :math:`-q\mu_{I^{-}}n_{I^{-}} \nabla\phi_{I^{-}} = \frac{3}{2}\left ( \frac{ - E_{red}(\mathbf{r_{c}})}{R_{L}} \right )`: split of the current between the ionic species;

*  :math:`-q\mu_{I^{-}_{3}} n_{I^{-}_{3}}\nabla\phi_{I^{-}_{3}} =  -\frac{1}{2}\left ( \frac{ - E_{red}(\mathbf{r}_{c})}{R_{L}} \right )`: split of the current between the ionic species;
*  :math:`\nabla\phi_{c} = 0`: no cationic current;


integral boundary conditions for conservation of ionic species:

*  :math:`\int_{\Omega} \left [ \frac{1}{3}n_{I^{-}}(\mathbf{r}) + n_{I^{-}_{3}}(\mathbf{r}) \right ] d\mathbf{r} = \left (\frac{1}{3}\bar{n}_{I^{-}} + \bar{n}_{I^{-}_{3}} \right )\Omega`: conservation of iodine ions within the cell;

*  :math:`\int_{\Omega} n_{c}(\mathbf{r}) d\mathbf{r}  =  \bar{n}_{c}\Omega`: conservation of cation within the cell;


where :math:`\Omega` is the volume of the cell, :math:`n_{\alpha}` the density
of charged species and the index :math:`\alpha` stands for cation (c),
iodide (I:math:`^{-}`), triiodide (I:math:`^{-}_{3}`) and electrons (e).
R:math:`_{L}` is the external load. The bias applied is equal to:

.. math::
   
   :label:
   
   \begin{equation}\label{pot1}
    V = \phi_{e} - E_{red}/q.
   \end{equation}

Ered is the redox potential. The redox potential can be evaluated using a Nernst approximation:

.. math::
   
   :label:
   
   \begin{equation}\label{redox11}
    E_{red} = E^{0}_{Pt} - \frac{kT}{2} ln \left ( \frac{n_{I^{-}_{3}}/n_{St}}{ (n_{I^{-}}/n_{St})^{3} } \right ).
   \end{equation}

Module DSC
----------------------

The DSC module is tagged as ``dssc``. In this part of the input file are set the name of the
simulation and the list of plotted variables:

::

  Module dssc 
    {
     name = dssc
     plot = (Potential, Density, Current, ContactCurrents)
     Solver linesearch
       {
        max_iterations = 300
        step_tolerance = 1e-4
        linear_solver
          {
           preconditioner = lu
          }
       }
     Physics
       {
        ...
       }
     Contact anode
       {
        ...
       }
     Contact cathode
       {
        ...
       }
    }

This section contains information about the ``Contacts``, parameters for the ``Solver`` and
the ``Physics`` sections.

Contacts
^^^^^^^^^^^^^^^^^^^^^

Information for the contacts is inserted in the ``Module`` section, in the subsections ``Contacts``. 
The two contacts are the photoanode and the cathode. The photoanode must be
the boundary of a region where :math:`\bf TiO_2` is present, on the contrary the cathode must be
the boundary of a region where the *electrolyte* is present.

::

  Contact anode
    {
     type = ohmic
     bias = $V[0.0]
    }
  Contact cathode
    {
     type = Pt
     load = $R[1e8]
    }

The anode is an ohmic contact, it contains only the sweep of the bias applied for the
electrochemical potential of the electrons. The cathode must contain the initial external
resistance (``load = R``). 

If a simulation of the cell under illumination is made the initial
value of ``load`` must be set to a high value, in the range of :math:`10^{6}-10^{8}\Omega cm^{2}`

in order to have
no current during the light intensity sweep. In case a simulation under dark conditions
is performed, where an external bias only is applied, ``load = 1.0``.

Physics
^^^^^^^^^^^^^^^^

For the generation:

|  ``generation = dssc_generation``

The ``Physics`` section must contain at least the setting of the generation module. The
set of parameters that can be defined for the entire device are enlisted in table :ref:`Physical parameters<dsc_parameters>`.

.. _dsc_parameters:

.. math::
   
   :label:
   
   \begin{table}[!h]
   \center
   \begin{tabular}{l|p{8cm}|l}
   \multicolumn{3}{c}{\textbf{Parameters}} \\
   \hline
    & \textbf{Poisson equation} & \\
   \hline
   \texttt{porosity} & Porosity of porous material &  \\
   \texttt{perm\_oxide} & Relative perm. of the TiO$_2$ &  \\
   \texttt{perm\_electrolyte} & Relative perm. of the electrolyte & \\
   \texttt{k\_dye} & Oxidized Dye regeneration rate constant & s$^{-1}$ \\
   \texttt{trap\_exp} & Trap exponential tail $a$& \\
   \texttt{trap\_DOS} & Trap effective density of states $N_t$& \\
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
   \texttt{rec\_non\_linearity} & Non-linearity exponent in the electron density recombination &  \\
   \end{tabular}
   \caption{Physical parameters that can be set for the model divided
   in subsets relative to different processes in the cell.}
   \label{table:dsc_parameters}
   \end{table}


| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 
| 



Generation module
^^^^^^^^^^^^^^^^^^^^^^^^^

The generation term is related to the flux of photons which reaches the active :math:`TiO_2`
regions and the dye present in the cell. We assume a simple Beer-Lambert exponential
decay for charge generation of the form:

..  math::
    
    :label:
    
    \begin{equation}\label{generation}
    G(\mathbf{r}) = \int^{\lambda_{2}}_{\lambda_{1}} \alpha(\lambda)
    \Phi (\lambda) e^{-\alpha(\lambda) \hat{n} \cdot \mathbf{r}}
    d\lambda
    \end{equation}

where :math:`\alpha` is the absorption coefficient (in :math:`\mu m^{-1}`
) of the chosen Dye, :math:`\Phi(\lambda)` the intensity of
the light power at wavelength :math:`\lambda` of the light source spectrum.
The part of the input file for the generation module:

:: 

  Module dssc_generation
    {
     regions = TiO2
     plot = (Distance, Generation)
     light_direction = (1, 0, 0)
     light_intensity = $x
     dye = N719
     
     Contact anode
       {
       }
    }

In the generation module must be specified:

* ``regions`` : the regions where we want that there is generation (where the Dye is present);
* ``plot if we want`` : to plot the generation;
* ``light_direction`` : the vector which fixes the direction from where the light comes;
* ``light_intensity`` : the light intensity;
* ``dye`` : the dye used in the cell;
* ``illumination_spectrum`` : the source spectrum, by default a 1.5 AM solar spectrum;

The light intensity is defined in a sweep. It can be set to reach 1 that means one Sun,
or a larger or smaller illumination intensity (0.1, 2.0, etc.).
There is another flag called ``illumination_spectrum`` that can be used if it is wanted
to change the spectrum profile :math:`\Phi` with another illumination source (for example if it is
assumed that the cell is under light concentration). The file used by default for F is the
standard 1.5 AM sun spectrum contained in the material database in the file **Sun1p5am** .

The last part of the generation is the definition of the Contact. This Contact tells to
the code which is the boundary region of the grid from where the light penetrates into
the device. The information given by the boundary and the light direction vector defines
the coming of light and direction of rays.

Device
----------------------

Device section for a DSC simulation:

::

  # Description of the device physical regions Device
  { 
   meshfile = <name_of_the_mesh>
   Region TiO2
     {
      material = TiO2mes
      porosity = 0.5
     }
   Region electrolyte
     {
      material = TiO2mes
      TiO2 = false
     }
  }

The Device section must contain the name of the mesh file ``meshfile = <name_of_the_mesh>`` .
For every region of the device it is specified the material file from the database (the
``TiO2mes`` contains standard parameters for both :math:`TiO_2` and electrolyte, so it can be used
for both regions). For every region must be specified if it contains ``TiO_2``, or electrolyte
or both. This can be done setting two flags called ``TiO2`` and electrolyte. 

If set *true* the material is present in the region, if set *false* it is not present. 

By default they are assumed both true (porous region). In the second region of the example showed here
there is electrolyte only and ``TiO2 = false`` is explicitly specified. In case both materials
are present (porous region) a porosity must be specified (in the range between 0, :math:`TiO_2`
only, and 1, electrolyte only). 

If one material is not present the porosity is automatically set to 0 or 1.

Sweep
-----------------

Three sweeps are needed for the plot of the entire I-V under illumination. The first sweep
is needed to make the transition from dark condition to full open-circuit condition under
illumination. The high value of the load **R** maintains the current to zero. Then a second
sweep is used to pass from open circuit condition to short circuit condition lowering the
external load from a high external load to a small one (``R = 1``). 
Finally, the voltage sweep compute the I-V characteristics under illumination. In case of dark simulation
(application of an external bias without illumination) the first and second sweep are not
needed. The external load for the cathode can be set directly to **R** = 1.

::

  Module sweep 
    {
     name = sweep_gen
     solve = (dssc_generation, dssc)
     variable = $x
     values = (0, 1e-9, 1e-8, 1e-7, 1e-6,
     1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1)
     plot_data = true
    }
  Module sweep 
    {
     name = sweep_R
     solve = dssc
     variable = $R
     values = ( 1000, 100, 10, 1 )
     plot_data = true
    }
  Module sweep 
    {
     name = sweep_V
     solve = dssc
     variable = $V
     values = (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.62, 0.64, 
     0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8)
     plot_data = true
    }
    
The intensity of illumination can be changed sweeping the value of **x** different to 1 (1 = 1 Sun of power).


Output
----------

The output that want to be plotted are enlisted in section Model within option. See 
:ref:`Table nodal<dsc_nodal>` , :ref:`Table elemental<dsc_elemental>` and :ref:`Table scalar<dsc_scalar>` .

..  _dsc_nodal:

..  math::
    
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



..  math::
    
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

..  _dsc_scalar:


    
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



.. _Solver_Dsc:

Solver
----------

Three sweeps are needed for the plot of the entire I-V under illumination. The first sweep
is needed to make the transition from dark conditions to full open-circuit conditions under
illumination. The high value of the load R maintains the current to zero. 
Then a second sweep is used to pass from open circuit condition to short circuit condition lowering the
external load from a high external load to a small one (R = 1). 
Finally, the voltage sweep compute the I-V characteristics under illumination. In case of dark simulation
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




..   rubric:: Footnotes

