.. _Dsc:


Simulation Dye Solar Cells
==============================

For a brief list of literature to understand Dye Solar cells (DSC) see [Kalyanasundaram]_ . For a review
of the model see [Gagliardi]_ .

The model consists in a set of drift-diflusion equations for the description of the
propagation of different ions and for the electrons coupled with Poisson equation:

..  math::
    :nowrap:
    :label:

    \begin{eqnarray}
    \nabla \cdot (\mu_{e}n_{e}\nabla\phi_{e}) & = & (G - R) \\
    \nabla \cdot (\mu_{I^{-}}n_{I^{-}}\nabla\phi_{I^{-}}) & = & -\frac{3}{2}(G - R) \\
    \nabla \cdot (\mu_{I^{-}_{3}}n_{I^{-}_{3}}\nabla\phi_{I^{-}_{3}}) & = & \frac{1}{2}(G - R) \\
    \nabla \cdot (\mu_{c}n_{c}\nabla\phi_{c}) & = & 0,
    \label{eq1}
    \end{eqnarray}

where :math:`\mu_\alpha` refers to carrier mobilities, :math:`n_\alpha` to concentrations and :math:`\phi_\alpha` to electrochemical
potentials. *R* is the recombination term and *G* the generation term due to illumination.
The Poisson equation to handle the internal potential drop:

..  math::
    :nowrap:
    :label:
    
    \begin{equation}\label{poisson}
    -\varepsilon\triangle \varphi =  q[n_{c} + N_{D}^{+} - n_{I^{-}} - n_{I_{3}^{-}} - (n_{e} -
    \bar{n}_{e})],
    \end{equation}

where :math:`N^+_D` is the amount of ionized dyes and it is equal to:

..  math::
    :nowrap:
    :label:

    \begin{equation}\label{dyeion}
    N_{D}^{+} = \frac{G}{k_{3}}
    \end{equation}
    
with *G* the generation term and :math`k_3` the rate constant of dye regeneration. The dielectric
constant, :math:`\epsilon` , of the mesoporous material is a mix the two dielectric functions of the semi-
conductor and the electrolyte. We use the Maxwell-Garnet model where the dielectric
function of the mixed medium becomes:

..  math::
    :nowrap:
    :label:

    \begin{equation}\label{diel}
    \varepsilon = \varepsilon_{s}\frac{\varepsilon_{e} + 2\varepsilon_{s} + 2\epsilon_{p}\varepsilon_{e} - 2\varepsilon_{s}\epsilon_{p}}
    {\varepsilon_{e} + 2\varepsilon_{s} -\epsilon_{p}\varepsilon_{e} +\varepsilon_{s}\epsilon_{p}}
    \end{equation}

where :math:`\epsilon_s` and :math:`\epsilon_e` are the dielectric constants of the semiconductor and the electrolyte,
respectively, and :math:`\epsilon_p` is the porosity of the medium. The recombination term depends
largely on the loss mechanisms at the electrolyte/oxide interface which follows the reac-
tion chain:

..  math::
    :nowrap:
    :label:

    \begin{eqnarray}
    I^{-} & \rightleftharpoons & I + e \\
    2I & \rightleftharpoons & I_{2} \\
    I_{2} + I^{-} & \rightleftharpoons & I^{-}_{3}.
    \label{reaction_loss}
    \end{eqnarray}

From the chemical path we can get a formula for the interface recombination (considering
that the flrst chemical reaction is the slow process):

..  math::
    :nowrap:
    :label:

    \begin{equation}\label{ricombinazione}
    R = k_{e} \left [  \left ( \frac{n_{e}}{\bar{n}_{e}} \right )^{\beta}\bar{n}_{e}\sqrt{\frac{n_{I^{-}_{3}}}{n_{I^{-}}}}
    - \bar{n}_{e}\sqrt{\frac{\bar{n}_{I^{-}_{3}}}{\bar{n}^{3}_{I^{-}}}} n_{I^{-}}\right
    ],
    \end{equation}  
    
where the electron rate :math:`k_e` is the recombination rate constant.

For the boundary conditions of the model we assume at the photoanode:
*  :math:`\phi_{e} = V` : electrochemical potential of electrons set with the voltage applied;
*  :math:`\nabla\phi_{I^{-}}  = 0` : no iodide current at the photoande;
*  :math:`\nabla\phi_{I_{3}^{-}} = 0` : no triiodide current at the photoande;
*  :math:`\nabla\phi_{c} = 0` : no cationic current;
*  :math:`\nabla\varphi = 0` : no charged layer at the photoanode;

at the cathode:

*  :math:`\nabla\phi_{e} = 0` : no electronic current at the cathode;
*  :math:`-q\mu_{I^{-}}n_{I^{-}} \nabla\phi_{I^{-}} = \frac{3}{2}\left ( \frac{ - E_{red}(\mathbf{r_{c}})}{R_{L}} \right )` : split of the current between the ionic species;
*  :math:`-q\mu_{I^{-}_{3}} n_{I^{-}_{3}}\nabla\phi_{I^{-}_{3}} =  -\frac{1}{2}\left ( \frac{ - E_{red}(\mathbf{r}_{c})}{R_{L}} \right )` : split of the current between the ionic species;
*  :math:`\nabla\phi_{c} = 0` : no cationic current;

integral boundary conditions for conservation of ionic species:

*  :math:`\int_{\Omega}$ $\left [ \frac{1}{3}n_{I^{-}}(\mathbf{r}) + n_{I^{-}_{3}}(\mathbf{r}) \right ] d\mathbf{r} = \left (\frac{1}{3}\bar{n}_{I^{-}} + \bar{n}_{I^{-}_{3}} \right )\Omega` : conservation of iodine ions within the cell;
*  :math:`\int_{\Omega} n_{c}(\mathbf{r}) d\mathbf{r}  =  \bar{n}_{c}\Omega` : conservation of cation within the cell;

where :math:`\omega` is the volume of the cell, :math:`n\alpha` the density of charged species and the index :math:`\alpha`
stands for cation (c), iodide ( :math:`I^{-}` ), triiodide ( :math:`I^{-}_3` ) and electrons (e). :math:`R_L` is the external
load. The bias applied is equal to:

..  math::
    :nowrap:
    :label:

    \begin{equation}\label{pot1}
    V = \phi_{e} - E_{red}.
    \end{equation}

:math:`E_{red}` is the redox potential. The redox potential can be evaluated using a Nernst approx-
imation:

..  math::
    :nowrap:
    :label:
    
    \begin{equation}\label{pot1}
    V = \phi_{e} - E_{red}.
    \end{equation}
