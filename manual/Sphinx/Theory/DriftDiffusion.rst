..   <marker>

.. _DriftDiffusionTheory:

Drift-diffusion simulation of electrons and holes
=================================================

Theory
----------------


The semi-classical transport simulation of electrons and holes is based on the drift-diffusion
approximation (see [Selberherr]_ ).

Beside the electric potential the electro-chemical potentials are used as variables such
that the system of PDEs to be solved reads as follows

.. math::
   :nowrap:
   :label:
   
   \begin{eqnarray}
   -\nabla(\varepsilon\nabla\varphi - \mathbf{P}) & = & -e(n - p - N_d^+ + N_a^-) \nonumber \\
   -\nabla(\mu_n n ( \nabla\phi_n + P_n \nabla T)  ) & = & R \\
   -\nabla(\mu_p p (\nabla\phi_p + P_p \nabla T) ) & = & -R \nonumber
   \end{eqnarray}
 
*P* is the electric polarization due to e.g. piezoelectric effects and *R* is the net 
recombination rate, i.e. recombination rate minus generation rate. :math:`P_n` and :math:`P_p` are the electron
and hole thermoelectric power, respectively. The models for the mobilities and the net
recombination rates can be specified in the **physical_model** sections as described in the
following.

..  index:: double:DriftDiffusion;Solution

Solution/Plot variables
-----------------------

The solution variables available for plotting and for interaction with other models are
given in :ref:`Table Plotting variables<dd_solutions>` .


Module options
---------------------

The following options influence the behaviour of the Drift-Diffusion module:

|  ``coupling = string`` defines which equations to solve. 
| 
|          The default is full, meaning that
           the full system consisting of the Poisson, electron continuity and hole continuity
           equations is solved. Other possible values are poisson (for equilibrium calculations),
           electrons or holes. For the last two cases local equilibrium is assumed such that
           :math:`\phi_n  =  \phi_p` .
 
|  ``enforce_local_charge_neutrality = bool`` 
| 
|        If set to true, local charge neutrality will
         be enforced by setting the charge density to zero. This may be useful for solving
         the Poisson equation involving only dielectrics.

|  ``guess_el_qfermi = double`` 
|  
|        If this option is set, then before resolving the system the
         given number will be set as a guess for the electron electro-chemical potential.

|  ``guess_hl_qfermi = double`` 
| 
|        If this option is set, then before resolving the system the
         given number will be set as a guess for the hole electro-chemical potential.

|  ``default_boundary_condition = string`` 
| 
|        With this option the user can control the default 
         boundary condition for the electric field on all external boundaries without explicit boundary model. 
         Possible values are zero field (default), or zero displacement.
         The two differ only in presence of electric polarization fields.

|  ``quadrature_rule = string`` 
| 
|        This option allows to chose between trapezoidal and Gauss
         type numeric integration rules. The default rule is gauss, but in some cases trapez
         may prevent density peaks near badly resolved material interfaces.

|  ``save_state = boolean`` 
|
|        If set to *true* the current solution will be written to a compressed 
         file after each solve. The file name follows the same rules as the result files,
         having file extension ``.tsv`` .

|  ``load_state = file`` 
| 
|       Reload a formerly saved solution. filename can be a filename (to
        reside in the current working directory), a relative or an absolute path to a file.
        The file needs to have been created with ``save_state = true`` .

|  ``solve_after_load = boolean`` 
| 
|       If set to true, the system will be solved after having
        reloaded a saved state. Otherwise it will not be solved, which is the default behaviour.

|warn|  
            Currently the reload of saved solutions *only works correctly when using the identical mesh* . 
            Otherwise there will be undefined behaviour or failure. 
            Future releases will relax this restriction.

Solver section
--------------------

The **Solver** section of the Drift-Diffusion module refers to a nonlinear solver .

Physics section
--------------------

The Physics block contains generic options for the bulk physical model and the definition
of submodels. The generic options are:

|  ``model = string`` 
| 
|    Specify the semiconductor model to be used. Possible values are default
     and simple. The former uses k.p to calculate the (strain corrected) band parameters.

|  ``thermal_simulation = string`` 
| 
|    If you are doing coupled electrothermal simulations,
     you have to specify the name of the thermal simulation providing the lattice temperature.

|  ``strain_simulation = string`` 
| 
|    If you are doing simulations on strained systems, you
     can specify the name of a strain simulation. The band parameters will then be
     calculated taking local strain corrections into account.

|  ``relax_polarization = double`` 
| 
|    With this option one can specify a relaxation factor
     for the electric polarization field. This can be useful if the amount of total electric
     polarization has to be treated as fitting parameter.

For the simple semiconductor model one has to provide conduction and valence band
edges and the effective density of states masses (or the effective density of states itself)
in the ``Region`` sections. The corresponding keywords are given in table 2.2 .
In the following we describe the submodels. All submodels can be restricted to a
subset of simulation regions.

Recombination/generation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the currently available generation/recombination models. Note
that all recombination models can be applied also to surfaces/interfaces.
Recombination models are controlled by means of recombination submodel blocks
inside the Physics section. Different recombination models having the same (or no)
options can be enabled in a single statement writing::

  recombination (model1, model2, ...) {}

**Shockley-Read-Hall (SRH) recombination**

The SRH recombination model can be enabled by defining a recombination submodel
of type ``srh`` .

**SRH** recombination is defined as follows:

..  math::
    :nowrap:
    :label:
    
    \begin{equation}
    R_{SRH} = \frac{np - n_i^2}{(n + n_i e^{E^*/k_BT})\tau_p +
    (p + n_i e^{-E^*/k_BT})\tau_n}
    \end{equation}

|

:math:`E^* = E_{trap} - (E_c + E_v)/2` is the trap level with respect to the midband energy. 

:math:`n_i` is the intrinsic carrier density, :math:`\tau_n` and :math:`\tau_p` are the recombination times. 

The parameters are taken from the material database. The recombination times are dependent on temperature and
doping density, e.g.

..  math::
    :nowrap:
    :label:
    
    \begin{eqnarray}
    \tau_n & = & \tau_n^0 \left(\frac{T}{T_0}\right)^{\alpha_n} e^{\beta(T/T0 - 1)} \\
    \tau_n^0 & = & \tau_{min,n} + \frac{\tau_{max,n} - 
    \tau_{min,n}}{1 + (N/N_{ref})^\gamma}
    \end{eqnarray}

where :math:`T_0` is the reference temperature (300 K). Table 2.3 shows the corresponding parameters 
for the material data files. The parameters for holes and electrons have to be
specified in an array, e.g. :math:`\tau_min = (1e-5, 3e-6)`

    

The recombination times and trap level can be overridden from the input file by using
the keywords of Table 2.4.



The SRH recombination model can be applied also to surfaces and interfaces. In this
case, you can provide the recombination velocities using the keywords ``rec_velocity_n``
and ``rec_velocity_p`` instead of ``tau_n`` and ``tau_p`` .

**Direct (radiative) recombination**

The direct recombination model can be enabled in the input file by by defining a
``recombination`` submodel of type ``direct`` .

Direct recombination is modeled as follows:

.. math::
   :nowrap:
   :label:

    \begin{equation}
    R_{direct} = C(np - n_i^2)
    \end{equation}

The material data file and the input file use the same keyword C for the parameter C. The
database value can be overridden from the input file as described for SRH recombination.

**Auger recombination**


The Auger recombination model can be enabled in the input file by defining a recombination
submodel of type ``auger`` .

Auger recombination is modeled by the following equation

.. math::
   :nowrap:
   :label:

    \begin{equation}
    R_{auger} = (C_nn + C_pp)(np - n_i^2)
    \end{equation}

with temperature dependent parameters

.. math::
   :nowrap:
   :label:

    \[
    C_{\{n,p\}}  =  \left(A + B\frac{T}{T_0} + C\left(\frac{T}{T_0}\right)^2\right)\left(1 + H e^{-\{n,p\}/N_0}\right)
    \]
 
The parameters A;B;C;H and :math:`N_0` are taken exclusively from the database. They are
different for :math:`C_n` and :math:`C_p` and have to be specified as arrays with keywords A, B, C, H, N0,
e.g. A = (1e-31, 1e-32). The calculated values for :math:`C_n` and :math:`C_p` can be overridden from
the input file by specifying values for the keywords :math:`C_n` and :math:`C_p` .


**Optical generation**

A very simple model for photoelectric generation of electron-hole pairs is implemented
in TIBERCAD. It is enabled by specifying a ``generation`` submodel of type optical The
model imposes a constant generation rate which has to be provided by the keyword G in
units of :math:`(cm*s)^{-1}` . 

|warn| 
          Note that the simulation usually should define a sweep on the value
          of G from 0 to the desired generation.

Thermoelectric power models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermoelectric power models are the same for electrons and holes. 

The keyword is  ``thermoelectric_power`` , i.e ::

  thermoelectric_power [type]
    {
    }
    
The model keyword can be ``constant`` (i.e. the thermoelectric powers are read from the
database) or ``diffusivity_model`` where the thermoelectric powers are computed by

.. math::
   :nowrap:
   :label:

    \begin{equation}
    P_n = - \frac{k_b}{q}\left( \frac{5}{2} + \frac{e \phi_n + E_c - e \varphi}{k_b T} \right)
    \end{equation}

.. math::
   :nowrap:
   :label:
    
    \begin{equation}
    P_p = \frac{k_b}{q}\left( \frac{5}{2} - \frac{e \phi_p + E_v - e \varphi}{k_b T} \right)
    \end{equation}

The default is :math:`P_n = P_p = 0`

Mobility models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The models to be used for electrons and holes can be defined in a single submodel
block or independently using two blocks. The corresponding keywords are mobility or
 ``electron_mobility`` and ``hole_mobility`` , i.e.

::

  mobility [type]
    {
    }
    
  electron_mobility [type]
    {
    }
    
  hole_mobility [type]
    {
    }
    
When using the first approach, both carriers will use the same model, and parameters 
provided in the input file will also be used by both carriers. When mixing the
different definitions, the blocks ``electron_mobility`` and ``hole_mobility`` will override
the common ``mobility`` block.

The default model is the constant mobility model. The parameters for the different
mobility models are needed for both electrons and holes. In the material files they are
specified with a common keyword in arrays, e.g.

..  math::
    :nowrap:
    :label:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l|c}
    \multicolumn{2}{c}{\textbf{Mobility Table}} \\
    \hline
    \textbf{mobility} & \textbf{constants}  \\
    \hline
    \hline
     & \texttt{electrons holes }  \\
    \hline
     &  \\
    \texttt{mu\_max} & (1400.0 , 250.0)  \\
    \texttt{exponent } & (1.0    , 2.1)  \\
    \hline
    \end{tabular}
    \caption{Mobility Model}
    \end{table}


Constant mobility model
-------------------------------

The constant mobility model (identifier ``constant`` ) assumes a mobility which depends
only on temperature by means of the following formula:

.. math::
   :nowrap:
   :label:

    \begin{equation}
    \mu_{const} = \mu_0 (T/T_0)^{-\gamma}
    \end{equation}

In the material data file :math:`\mu_0` and :math:`\gamma`  
have to be specified with the keywords ``mu_max`` and ``exponent``. 
:math:`\mu_0` can be overridden from the ``physical_model`` section using the keyword
mu or from the Region sections using the keywords ``mu_e`` and ``mu_h`` .

Doping dependent mobility model
-------------------------------------

The doping dependent mobility model (identifier ``doping_dependent`` ) implements two
models for mobility depending on the total doping density and the temperature. The
model that is used depends on the value of the ``mobility_formula`` parameter.

Model by Masetti et al. [4]

The model by Masetti et al. is identified by ``mobility_formula`` = 1. It uses the following
formula:

.. math::
   :nowrap:
   :label:
    

    \begin{align}
    \mu = \mu_{min,1} * \mathrm{e}^{-P_c/N} +
    \frac{\mu_{const} - \mu_{min,2}}{1 + (N/C_r)^\alpha} -
    \frac{\mu_1}{1 + (C_s/N)^\beta}
    \end{align}


where N is the total doping density and :math:`\mu_{const}` the mobility obtained from the constant
mobility model. The parameters are specified in the material file as given in Table 2.5.



Model by Arora [5]
The model by Arora is identified by ``mobility_formula`` = 2. It reads:

.. math::
   :nowrap:
   :label:

    \begin{align}
    \mu = \mu_{min} + & \frac{\mu_d}{1+(N/N_0)^{A^*}} \\
    \end{align}
    
with

.. math::
   :nowrap:
   :label:

    \begin{align}
    \mu_{min} = A_{min}(T/T_0)^{\alpha_m}, & \quad \mu_d = A_d(T/T_0)^{\alpha_d} \nonumber \\
    N_0 = A_N(T/T_0)^{\alpha_N}, & \quad A^* = A_a(T/T_0)^{\alpha_a} \nonumber 
    \end{align}

The parameters are given in table at the end of the Chapter.

|

**Field dependent mobility model**

The field dependent mobility model describes the degradation of mobility at high driving
fields. It is identified by the identifier field_dependent. The electric field component
in direction of the current 
ow or the gradient of the electro-chemical potential can be
chosen as driving force:

  ``driving_force = efield | grad_fermi | field_parameter``

The default driving force is the gradient of the corresponding electro-chemical potential :math:`\nabla\phi` .
``field_parameter`` uses a field parameter given by :math:`\sqrt{E\cdot\nabla\phi}`  as driving force.

The model is based on the Caughey-Thomas model, refined by Canali [6]:

.. math::
   :nowrap:
   :label:

    \begin{equation}
    \mu = \frac{\mu_{lowfield}}{\left(1 + \left(\frac{\mu_{lowfield} |\mathbf{E}|}{v_{sat}}\right)^\beta \right)^{1/\beta}}
    \end{equation}
    
with

.. math::
   :nowrap:
   :label:

    \[
    \beta = \beta_0(T/T_0)^b 
    \]

:math:`|E|` is the modulus of the driving field, :math:`\mu_{lowfield}` is the low-field mobility. For the latter
one can specify the model to be used using the parameter ``lowfield_model`` . As default
the doping dependent model is used.
There are two models for vsat, identified with ``Vsat_Formula = 1`` and 2. 

Formula 1 reads

.. math::
   :nowrap:
   :label:

    \[
    v_{sat} = v_{sat,0} (T/T_0)^{-\gamma}
    \]
 

Formula 2 reads

.. math::
   :nowrap:
   :label:

    \[
    v_{sat} = \max(A_{vsat} - B_{vsat} (T/T_0), v_{min})
    \]

The parameters for the field dependent mobility model are summarized in Table 2.7.



Polarization models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For simulations involving materials with nonzero electric polarization (such as nitrides)
it is important to include the effect of polarization. This is done by specifying the models
for spontaneous (pyro-) and piezoelectric polarization using the keywords polarization
with the types ``pyro`` and ``piezo`` ::

  polarization pyro {}
  polarization piezo {}

As for all models, if they do not have individual options, they can be specified together
by writing **polarization (pyro, piezo) {}** .

|

**Spontaneous (pyro-) polarization**

The spontaneous polarization model imposes a constant electric polarization P along
the symmetry-braking direction of the crystal. Crystals with wurtzite structure like Nitrides
have strong piezoelectric fields along the c-direction. 
The value of the polarization usually is taken from the database, but it can be overridden from the input file by specifying
the option Pz, meaning the value of the spontaneous polarization along c-direction ([0001]). 
Alternatively, one can specify explicitly a polarization vector using the option
**P = (Px, Py, Pz)** . This is useful to impose an arbitrary constant polarization field.

|

**Piezopolarization**

The piezoelectric polarization is strain induced and given by

.. math::
   :nowrap:
   :label:

    \begin{equation}
    P^{pz} = e_{ikl}\varepsilon_{kl}
    \end{equation}

where :math:`\varepsilon_{kl}` is the strain tensor. The piezoelectric moduli :math:`e_{ikl}` are stored in the database.
The strain is obtained from the simulation specified in the ``Physics`` section, but it can
be overridden by providing a name for the strain simulation inside the polarization block
using the ``strain_simulation`` option.

particle density
^^^^^^^^^^^^^^^^^^^^^^^^^^

Details for the calculation of the electron and hole densities can be given in the particle_density
submodel. Its options are:

|  ``particle = string``   The particle this model is describing. Can be ``electron`` or ``hole`` .

|  ``statistics = string``   The statistics. Can be ``fermidirac`` (default) or ``boltzmann`` .

|  ``quantum_density = string`` The name of a quantum density simulation. 
| 
|       This will use the quantum mechanical particle density in the regions it was calculated.

If a quantum density is used, then it is useful to define also an embracing region
where the model gradually switches from a fully classical to a fully quantum density.
The options for the embracing are specified in a block with keyword ``embracing`` . It
accepts the following options:

|  ``embracing_length = double`` 
| 
|      When the domain of the quantum simulation is smaller
       than the domain of the full simulation, the boundary conditions for the Schroedinger
       equation will disturb the transfer from classical to quantum density. By defining an
       embracing region of a certain extension (specified in meters), a gradual transition
       from classical to quantum density will be done instead of an abrupt one, using as
       effective density :math:`\rho = x\cdot\rho_{\mathrm{quantum}} + (1-x)\cdot\rho_{\mathrm{classical}}` . 
       The default is no embracing region at all (zero extension).

|  ``cutoff = double`` 
| 
|      If an embracing region is used, a part of this region near the boundary
       of the quantum region can be cut off so that only the classical density is considered
       in that part. ``cutoff`` is specified as a percentage of the embracing length and should
       therefore be between 0.0 and 1.0.

|  ``plot_embracing_region = bool`` 
| 
|      Whereas the automatic creation of the embracing region 
       in 1D is a very simple task, it is a more difficult one in higher dimensions. 
       By setting this flag to true, the embracing region and the mixing coefficient x will be
       plotted for a visual control of the quality of the embracing region. 
       The default is ``false`` .

Trap models
^^^^^^^^^^^^^^^^^^^^^^^^

Currently single level traps are implemented in TiberCAD. Traps can be normally neutral
or normally charged electron or hole traps, or a fixed charge. Common options for all
models are

|  ``type = string`` (If not provided as second keyword). 

|       The type of traps. One of eNeutral, hNeutral, donor, acceptor or fixed_charge.

|  ``Nt = double`` The trap density in cm 3 (or cm 2 for surface traps).

|  ``Et = double`` The trap level with respect to a reference energy.

|  ``reference = string`` The reference energy. The default is m. 

The trap energy in this case is given as :math:`E_{trap} = E_{midgap} + Et` . 

Other possible values are :math:`E_{trap} = E_c - Et` or :math:`E_{trap} = E_v + Et` .

|   ``eNeutral`` The trapped electron density is given by

..  math::
    :nowrap:
    :label:
    
    \[
    n_t = \frac{N_t}{1 + \exp(\frac{E_{trap} - E_{F,n}}{k_BT})}
    \]

|   ``hNeutral`` The trapped hole density is given by

..  math::
    :nowrap:
    :label:
    
    \[
    p_t = \frac{N_t}{1 + \exp(-\frac{E_{trap} - E_{F,p}}{k_BT})}
    \]

|   ``donor`` The density of ionized traps is given by

..  math::
    :nowrap:
    :label:

    \[
    N^+_t = N_t - \frac{N_t}{1 + \exp(\frac{E_{trap} - E_{F,n}}{k_BT})}
    \]

|   ``acceptor`` The density of ionized traps is given by

..  math::
    :nowrap:
    :label:
   
    \begin{equation}
    N^-_t = N_t - \frac{N_t}{1 + \exp(-\frac{E_{trap} - E_{F,p}}{k_BT})}
    \end{equation}

If traps are specified, the total charge density in the Poisson equation is modified to
include the charged trap densities:

.. math::
   :nowrap:
   :label:

   \begin{equation}
   \rho = e\left(p - n + N^+_D - N^-_A - \sum n_t + \sum p_t + \sum N^+_t - \sum N^-_t \right)
   \end{equation}

Additionally, each trap induces a SRH recombination term of the form

.. math::
   :nowrap:
   :label:

   \begin{equation}
   R_t = N_t \frac{v_{th}^n\sigma^nv_{th}^p\sigma^p(np -n_i^2)}{v_{th}^n\sigma^n(n+n_1) + v_{th}^p\sigma^p(p+p1)}
   \end{equation}

where :math:`\sigma^{n,p}` are the capture cross sections, :math:`v_{th}^{n,p}` the thermal velocities and (for Boltzmann statistics)

.. math::
   :nowrap:
   :label:

   \begin{equation}
   n_1 = n_{i,eff} \exp(E_{trap}/k_BT),\quad p_1 = n_{i,eff} \exp(-E_{trap}/k_BT)
   \end{equation}

Boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Boundary conditions are implemented for ohmic contacts, Schottky contacts, free surfaces
and interfaces. Contacts are boundary models that allow a nonzero normal electrical
current. The applied voltage is specified with the option ``voltage`` . 

A variable can be assigned to this, using the $-syntax. On ohmic or schottky contacts one can define surface
recombination velocities for electrons and holes using the options ``rec_velocity_e`` and 
``rec_velocity_h`` . This will impose Robin type boundary conditions for the continuity
equations of the form

.. math::
   :nowrap:
   :label:

    \begin{eqnarray}
    -\nabla(\mu_n n \nabla\phi_n ) & = & v_n(n - n_0) \\
    -\nabla(\mu_p p (\nabla\phi_p + P_p \nabla T) ) & = & -v_p(p - p0)
    \end{eqnarray}

The options ``zero_field`` , ``zero_grad_fermi_e`` and ``zero_grad_fermi_h`` can be used,
which when set to ``true`` will impose zero normal electric field and zero normal gradient
of the electron and hole electro-chemical potential, respectively. The latter are special
cases of surface recombination velocities ( :math:`v_{r}ec = 0` ).

Contacts are defined by blocks with keyword Contact, for example::

  Contact anode 
    {
     type = ohmic
     [regions = (anode1, anode2)]
     voltage = $Vd
    }
    
An area factor can be specified for contacts using the keyword ``area_factor`` . The
contact current will be multiplied by this factor.

For interfaces and surfaces, the same syntax can be used (optionally one can use
the keywords ``Interface`` or ``Boundary`` ), however they do usually not need to be defined
explicitly.

**Ohmic contact**

The ohmic contact (identifier ohmic) has no further parameters.

**Schottky contact**

A Schottky contact (identifier ``schottky`` ) has the additional parameter ``barrier`` , which
signifies the energy difference between the semiconductor band edge and the fermi energy
in the metal. As default, the barrier is taken with respect to the conduction band. By
specifying ``band = v`` the barrier can be imposed with respect to the valence band (p-
type contact). Alternatively, the metal work function can be defined using the keyword
work_function. Note, however, that its value has to be aligned with the band energies
given in the material files for the other materials. The ``fixed_barrier`` controls the
behaviour of the barrier height for strained materials. If it is set to true, the barrier will
be independent of strain (default behaviour). If it is set to ``false`` , the given barrier is
used as barrier for the unstrained case and will depend on strain during simulation. If
the metal work function is specified, the barrier will be strain dependent as default.

  |warn| 
            if the contact is touching different materials, one should specify the work
            function instead of the barrier.

Thermionic emission is by default switched on, but can be disabled by specifying
``thermionic_emission = false`` .

|

**Interface/surface model**

The free surface or interface model (identifier interface) can include surface charges
due to traps and surface recombination. Their definition can be found in section (see above).

Each trap model will induce automatically a SRH recombination model as in the bulk
case.

Schroedinger/Poisson/Drift-Diffusion calculations
-----------------------------------------------------

TIBERCAD is able to do selfconsistent Schroedinger-Poisson or Schroedinger-Drift-Diffusion
calculations. For this purpose, quantum_density has to be specified for at least one of
the carriers, and a selfconsistent simulation should be defined in the Selfconsistent
block. The following options { to be specified in the Physics section { control the
behaviour of the selfconsistent simulation.

|  ``use_density_predictor = bool`` 

When set to true, a predictor-corrector scheme will
be adopted in the selfconsistent cycle. The Poisson/Drift-Diffusion solver does not
just take the particle densities as given by the Schroedinger calculation, but it will
assume a dependency of the density on the potentials of the form

.. math::
   :nowrap:
   :label:

    \begin{equation}
    \rho(\varphi, \phi_n, \phi_p) = \frac{\rho_{\mathrm{quantum}}(\varphi^0, \phi_n^0, \phi_p^0)}{\rho_{\mathrm{classical}}(\varphi^0, \phi_n^0, \phi_p^0)}\rho_{\mathrm{classical}}(\varphi, \phi_n, \phi_p)
    \end{equation}

|

where :math:`(\varphi^0, \phi_n^0, \phi_p^0)` are the potentials for which the quantum density was calculated. use_density_predictor = true is the preferred method for selfconsistent
Schroedinger-Poisson/Drift-Diffusion calculations and is enabled by default.

Example
--------------

The following example shows the Drift-Diffusion module definition for a pn junction.

----

::

  Module driftdiffusion
    {
     name = dd
     #regions = (pside, nside)
     Solver linesearch 
       {
       }
     Physics
       {
        recombination srh {}
        
        mobility doping_dependent {}
        
        Contact anode 
          { 
           voltage = $Vd 
          }

          Contact cathode 
            { 
             voltage = 0 
            }
       }
    }

----

Listing 3: Models section for drift-diffusion

..  _dd_solutions :

..  math::
    :nowrap:
    :label:
    
     \begin{table}[!ht]
     \center
     \begin{tabular}{l|c|l}
     \multicolumn{3}{c}{\textbf{Solution Table}} \\
     \hline
     \textbf{Keyword}  & \textbf{Description} & \textbf{Units}  \\
     \hline
     \hline
     \texttt{Ec} & Conduction band edge & eV \\
     \texttt{Ev} & Valence band edge & eV \\
     \texttt{eQFermi} & Electro-chemical potential of electrons & eV ($-e\phi_n$) \\
     \texttt{hQFermi} & Electro-chemical potential of holes & eV ($-e\phi_p$) \\
     \texttt{Ec0} & Conduction band edge without electric potential & eV \\
     \texttt{Ev0} & Valence band edge without electric potential & eV \\
     \texttt{Eg} & Band gap & eV \\
     \texttt{ConductionBands} & Minimal of all conduction bands & eV \\
     \texttt{ValenceBands} & Maximum of all valence bands & eV \\
     \texttt{ElPotential} & Electric potential & V \\
     \texttt{eDensity} & Electron density & cm$^{-3}$ \\
     \texttt{hDensity} & Hole density & cm$^{-3}$ \\
     \texttt{eMobility} & Electron mobility & cm$^2V^{-1}s^{-1}$ \\
     \texttt{hMobility} & Hole mobility & cm$^2V^{-1}s^{-1}$ \\
     \texttt{eConductivity} & Electron conductivity & S/cm \\
     \texttt{hConductivity} & Hole conductivity & S/cm \\
     \texttt{ElField} & Electric Field & Vcm$^{-1}$ \\
     \texttt{Polarization} & Electric Polarization & Cm$^{-2}$ \\
     \texttt{CurrentDensity} & Total current density & Acm$^{-2}$ \\
     \texttt{eCurrentDensity} & Electron current density & Acm$^{-2}$ \\
     \texttt{hCurrentDensity} & Hole current density & Acm$^{-2}$ \\
     \texttt{IonizedDonors} & Ionized donor density & cm$^{-3}$ \\
     \texttt{IonizedAcceptors} & Ionized acceptor density & cm$^{-3}$ \\
     \texttt{IonizedElectronTraps} & Ionized electron trap density & cm$^{-3}$ \\
     \texttt{IonizedHoleTraps} & Ionized hole trap density & cm$^{-3}$ \\
     %\texttt{charge\_density} & Total charge density & cm$^{-3}$ \\
     \texttt{eThElPower} & Electron thermoelectric power & V K$^{-1}$ \\
     \texttt{hThElPower} & Hole thermoelectric power & V K$^{-1}$ \\
     \texttt{NetRecombination} & %\begin{minipage}[t]{8cm}
     %
     %The net recombination rate for each recombination model and the total rate
     %\end{minipage}& cm$^{-3}$s$^{-1}$
     The net recombination rate for each recombination/generation \\
       &  model and the total rate    cm$^{-3}$s$^{-1}$ \\
     \texttt{ContactCurrent} & The electric current on each contact & A/cm$^{3-d}$
     \end{tabular}
     \caption{Solution/Plot variables}
     \label{table:dd_solutions}
     \end{table}

|

..  math::
    :nowrap:
    :label:
    
     \begin{table}[!ht]
     \center
     \begin{tabular}{l||l}
     \multicolumn{2}{c}{\textbf{Semiconductor Table}} \\
     \hline
     \textit{keyword} & \textit{description} \\
     \hline \hline
     \texttt{Ec} & conduction band edge (eV) \\
     \texttt{Ev} & valence band edge (eV) \\
     \texttt{m\_dos\_e} & conduction band effective DOS mass (m$_e$) \\
     \texttt{m\_dos\_h} & valence band effective DOS mass (m$_e$) \\
     \texttt{Nc} & conduction band effective DOS (cm$^-3$) \\
     \texttt{Nv} & valence band effective DOS (cm$^-3$) \\
     \end{tabular}
     \caption{Parameters for the simple semiconductor model}
     \label{table:simple_sc}
     \end{table}

|

..  math::
    :nowrap:
    :label:
    
     \begin{table}[!ht]
     \center
     \begin{tabular}{l||l|l}
     \multicolumn{2}{c}{\textbf{SRH Table}} \\
     \hline
     \textit{parameter} & \textit{keyword} \\
     \hline\hline
     $\tau_{min}$ & \texttt{taumin} \\
     $\tau_{max}$ & \texttt{taumax} \\
     $N_{ref}$ & \texttt{Nref} \\
     $\gamma$ & \texttt{gamma} \\
     $E^*$ & \texttt{Etrap} \\
     $\alpha$ & \texttt{Talpha} \\
     $\beta$ & \texttt{Tcoeff} \\
     \end{tabular}
     \caption{SRH material data file parameters}
     \label{table:srh_params_db}
     \end{table}

|

..  math::
    :nowrap:
    :label:
    
     \begin{table}[!ht]
     \center
     \begin{tabular}{l||l}
     \multicolumn{2}{c}{\textbf{SRH parameters Table}} \\
     \hline
     \hline
     $tau_n$ & \texttt{tau\_n} \\
     $tau_p$ & \texttt{tau\_p} \\
     $E^*$ & \texttt{E\_t}
     \end{tabular}
     \caption{SRH input file parameters}
     \label{table:srh_params_input}
     \end{table}
    
|

.. math::
   :nowrap:
   :label:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l||l|l}
    \multicolumn{2}{c}{\textbf{Mobility Model Table}} \\
    \hline
    \textit{parameter} & \textit{keyword}  \\
    \hline\hline
    $\mu_{min,1}$ &  \verb+mumin1+ \\
    $\mu_{min,2}$ & \verb+mumin2+  \\
    $\mu_1$ & \verb+mu1+  \\
    $P_c$ & \verb+Pc+ \\
    $C_r$ & \verb+Cr+ \\
    $C_s$ & \verb+Cs+ \\
    $\alpha$ & \verb+alpha+ \\
    $\beta$ & \verb+beta+ \\
    \end{tabular}
    \caption{Data file parameters for the mobility model by Masetti et al.}
    \label{table:mobility_masetti}
    \end{table}
    
|

.. math::
   :nowrap:
   :label:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l||l|l}
    \multicolumn{2}{c}{\textbf{Arora Model Table}} \\
    \hline
    \textit{parameter} & \textit{keyword} \\
    \hline\hline
    $A_{min}$ &  \verb+mumin+ \\
    $A_d$ &  \verb+mud+ \\
    $A_N$ &  \verb+N0+ \\
    $A_a$ &  \verb+A+ \\
    $\alpha_m$ &  \verb+am+ \\
    $\alpha_d$ &  \verb+ad+ \\
    $\alpha_N$ &  \verb+aN+ \\
    $\alpha_a$ &  \verb+aA+ \\
    \end{tabular}
    \caption{Data file parameters for the mobility model by Arora.}
    \label{table:mobility_arora}
    \end{table}
    
|

.. math::
   :nowrap:
   :label:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l||l|l}
    \multicolumn{3}{c}{\textbf{Mobility Dependence Table}} \\
    \hline
    \textit{parameter} & \textit{keyword} \\
    \hline\hline
    $\beta_0$ &  \verb+beta0+ \\
    $b$ &  \verb+betaexp+ \\
    $v_{sat,0}$ &  \verb+vsat0+ \\
    $\gamma$ &  \verb+vsatexp+ \\
    $A_{vsat}$ &  \verb+A_vsat+ \\
    $B_{vsat}$ &  \verb+B_vsat+ \\
    $v_{min}$ &  \verb+vsat_min+ \\
    \end{tabular}
    \caption{Data file parameters for the mobility model by Arora.}
    \label{table:mobility_field_dep}
    \end{table}


.. rubric:: Footnotes

.. [#] the Default value is given in brackets.
.. [#] the linear tolerance gets automatically decreased after each nonlinear step.


..   </marker>

    
.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%



