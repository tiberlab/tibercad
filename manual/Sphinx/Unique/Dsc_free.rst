


Simulation Dye Solar Cells
========================================


Introduction
------------------------------

For a brief list of literature to understand Dye Solar cells (DSC) see [Kalyanasundaram]_ . 

For a review of the model see [Gagliardi]_ .

The model consists in a set of drift-diffusion equations for the propagation of ions
and electrons coupled with Poisson equation:

.. math::
   :label: dsc_eq_system

    \nabla \cdot (\mu_{e}n_{e}\nabla\phi_{e}) & =  (G - R) \\
    \nabla \cdot (\mu_{I^{-}}n_{I^{-}}\nabla\phi_{I^{-}}) & =  -\frac{3}{2}(G - R) \\
    \nabla \cdot (\mu_{I^{-}_{3}}n_{I^{-}_{3}}\nabla\phi_{I^{-}_{3}}) & =  \frac{1}{2}(G - R) \\
    \nabla \cdot (\mu_{c}n_{c}\nabla\phi_{c}) & =  0,

where :math:`\mu_{\alpha}` refers to carrier mobilities, :math:`n_{\alpha}` to charge concentrations and :math:`\phi_{\alpha}` to 
electrochemical potentials. **R** is the recombination term and **G** the generation term due to illumination. 
In order to take into account the trap density we use a density dependent
mobility developed from multi-trapping model:

.. math::
   :label: dsc_eq_difftraps

    \mu_{e}(n_{e}) = \mu_{0} \left ( \frac{n_{e}}{N_{t}} \right )^{\frac{1 - a}{a}}


where :math:`a` is the trap exponent, :math:`N_t` the trap density and :math:`\mu_0` a constant. The energy trap
density is assumed to form an exponential tail below the conduction band of the semiconductor:

.. math::
   :label: dsc_eq_denstrap
   
    g_{T}(E) = \frac{a N_t}{kT} e^{\frac{-a E}{kT}}.


The Poisson equation to handle the internal potential drop:

.. math::
   :label: dsc_eq_poisson
   
    -\varepsilon\triangle \varphi =  q[n_{c} + N_{D}^{+} - n_{I^{-}} - n_{I_{3}^{-}} - (n_{e} - \bar{n}_{e})],


where N+D is the amount of ionized dyes and it is equal to:

.. math::
   :label: dsc_eq_dyeion
   
    N_{D}^{+} = \frac{G}{k_{3}}


with **G** the generation term and :math:`k_{3}` the rate constant of dye regeneration. The dielectric
constant, :math:`\varepsilon` , of the mesoporous material is a mix the two dielectric functions of the semiconductor and the electrolyte. 
We use the Maxwell-Garnet model where the dielectric
function of the mixed medium becomes:

.. math::
   :label: dsc_eq_diel
   
    \varepsilon = \varepsilon_{s}\frac{\varepsilon_{e} + 2\varepsilon_{s} + 2\epsilon_{p}\varepsilon_{e} - 2\varepsilon_{s}\epsilon_{p}}
    {\varepsilon_{e} + 2\varepsilon_{s} -\epsilon_{p}\varepsilon_{e} +\varepsilon_{s}\epsilon_{p}}


where :math:`\varepsilon_{s}` and :math:`\varepsilon_{e}` are the dielectric constants of the semiconductor and the electrolyte,
respectively, and :math:`\epsilon_{p}` is the porosity of the medium. The recombination term depends
largely on the loss mechanisms at the electrolyte/oxide interface which follows the reaction chain:

.. math::
   :label: dsc_eq_reaction_loss
   
   I^{-} & \rightleftharpoons  I + e \\
   2I & \rightleftharpoons  I_{2} \\
   I_{2} + I^{-} & \rightleftharpoons  I^{-}_{3}.


From the chemical path we can get a formula for the interface recombination (considering
that the first chemical reaction is the slow process):

.. math::
   :label: dsc_eq_ricombinazione
   
    R = k_{e} \left [  \left ( \frac{n_{e}}{\bar{n}_{e}} \right )^{\beta}\bar{n}_{e}\sqrt{\frac{n_{I^{-}_{3}}}{n_{I^{-}}}}
    - \bar{n}_{e}\sqrt{\frac{\bar{n}_{I^{-}_{3}}}{\bar{n}^{3}_{I^{-}}}} n_{I^{-}}\right],


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
iodide (I\ :math:`^{-}`), triiodide (I\ :math:`^{-}_{3}`) and electrons (e).
R\ :math:`_{L}` is the external load. The bias applied is equal to:

.. math::
   :label: dsc_eq_pot1
   
    V = \phi_{e} - E_{red}/q.

Ered is the redox potential. The redox potential can be evaluated using a Nernst approximation:

.. math::
   :label: dsc_eq_redox11
   
    E_{red} = E^{0}_{Pt} - \frac{kT}{2} ln \left ( \frac{n_{I^{-}_{3}}/n_{St}}{ (n_{I^{-}}/n_{St})^{3} } \right ).





..   rubric:: Footnotes


.. |TiO2| replace:: :math:`{\rm TiO_2}`
