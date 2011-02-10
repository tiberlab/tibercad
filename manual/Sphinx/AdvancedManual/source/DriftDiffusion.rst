.. _DriftDiffusion:

Drift-diffusion simulation of electrons and holes
=================================================

Theory
------

The semi-classical transport simulation of electrons and holes is based on the drift-
diffusion approximation (see [Cit-C]_ ).

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

..  index:: single:Solution

Solution/Plot variables
-----------------------

The solution variables available for plotting and for interaction with other models are
given in :ref:`Table61` .

..  index:: single:Models

Models section
--------------

The **Models** section looks as given in :ref:`Listing3`
The **options** block supports three special keywords:

* **save_state = boolean** If set to **true** the current solution will be written to a compressed 
file after each solve. The file name follows the same rules as the result files,
having file extension **.tsv**.

* **load_state = file** Reload a formerly saved solution. *filename* can be a filename (to
reside in the current working directory), a relative or an absolute path to a file.
The file needs to have been created with **save_state = true**.

* **solve_after_load = boolean** If set to **true**, the system will be solved after having
reloaded a saved state. Otherwise it will not be solved, which is the default behaviour.

It is important that currently the reload of saved solutions *only works correctly when
using the identical mesh*. Otherwise there will be undefined behaviour or failure. Future
releases will relax this restriction.

The **physical_model** sections can be omitted. In this case default models are used,
namely no recombination/generation for the recombination models and constant mobility
for the mobility models. There can be more than one recombination model.

.. _srh:

..  index:: pair:Recombination;srh

Recombination/generation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the currently available generation/recombination models. Note
that all recombination models can be applied also to surfaces/interfaces.

**Shockley-Read-Hall (SRH) recombination**

The SRH recombination model can be enabled in the input file by::

  physical_model recombination

    {
     type = srh
     ...
    }

SRH recombination is defined as follows:

.. math::
   :nowrap:
   :label:
   
   \begin{equation}
   R_{SRH} = \frac{np - n_i^2}{(n + n_i e^{E^*/k_BT})\tau_p +
   (p + n_i e^{-E^*/k_BT})\tau_n}
   \end{equation}

----
   
.. _Listing3: 

Listing 3

::

  model driftdiffusion
    {
     options
       {
        simulation_name = whatever_you_want
        physical_regions = (2, 3, 4)
       }
     physical_model recombination
       {
        type = model_to_be_used
        ...
       }
     physical_model generation
       {
        type = model_to_be_used
        ...
       }
     physical_model electron_mobility
       {
        type = model_to_be_used
       }
     physical_model hole_mobility
       {
        type = model_to_be_used
       }
     physical_model thermoelectric_power
       {
        type = model_to_be_used
       }
    }
  
Listing 3: Models section for drift-diffusion

----

:math:`E* = E_{trap}-(E_c + E_v)/2` is the trap level with respect to the midband energy. ni is the
intrinsic carrier density, :math:`\tau_n` and :math:`\tau_p` are the recombination times. The parameters are taken
from the material database. The recombination times are dependent on temperature and
doping density, e.g.

.. math::
   :nowrap:
   :label:
   
   \begin{eqnarray}
   \tau_n & = & \tau_n^0 \left(\frac{T}{T_0}\right)^{\alpha_n} e^{\beta(T/T0 - 1)} \\
   \tau_n^0 & = & \tau_{min,n} + \frac{\tau_{max,n} - 
   \tau_{min,n}}{1 + (N/N_{ref})^\gamma}
   \end{eqnarray}
   
where :math:`T_0` is the reference temperature (300 K). :ref:`Table62` shows the corresponding pa
rameters for the material data files. The parameters for holes and electrons have to be
specified in an array, e.g. ``taumin = (1e-5, 3e-6)``

The recombination times and trap level can be overridden from the input file by using
the keywords of :ref:`Table63` in the appropriate **physical_model** section or in the **Region**
section (the latter overrides the former).

The SRH recombination model can be applied also to surfaces and interfaces. In this
case, you can provide the recombination velocities using the keywords **rec_velocity_n**
and **rec_velocity_p** instead of **tau_n** and **tau_p**.

.. index:: pair:recombination;radiative

**Direct (radiative) recombination**

The direct recombination model can be enabled in the input file by::

  physical_model recombination
    {
     type = direct
     ...
    }
    
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

The Auger recombination model can be enabled in the input file by6.3. MODELS SECTION 65
physical_model recombination::

  {
   type = auger
   ...
  }
  
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

The parameters *A,B,C,H* and :math:`N_0` are taken exclusively from the database. They are
different for :math:`C_n` and :math:`C_p` and have to be specified as arrays with keywords A, B, C, H, N0,
e.g. A = (1e-31, 1e-32). The calculated values for :math:`C_n` and :math:`C_p` can be overridden from
the input file by specifying values for the keywords Cn and Cp.

..  index:: pair:generation;optical

**Optical generation**

A very simple model for photoelectric generation of electron-hole pairs is implemented
in TIBERCAD. It is enabled by specifying::

  physical_model generation
    {
     type = optical
     G = ...
    }

The model imnposes a constant generation rate which has to be provided by the keyword
G in units of :math:`(cm*s)^{-1}`

Thermoelectric power models
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermoelectric power models are the same for electrons and holes. The keyword is
**thermoelectric_power** ::

  physical_model thermoelectric_power
    {
     model = ...
    }
  
The model keyword can be **constant** (i.e. the thermoelectric powers are read from the
database) or **diffusivity_model** where the thermoelectric powers are computed by

.. math::
   :nowrap:
   :label:
   
   \begin{equation}
   P_n = - \frac{k_b}{q}\left( \frac{5}{2} + \frac{e \phi_n + E_c - e \varphi}{k_b T} \right)
   \end{equation}

   \begin{equation}
   P_p = \frac{k_b}{q}\left( \frac{5}{2} - \frac{e \phi_p + E_v - e \varphi}{k_b T} \right)
   \end{equation}

The default is :math:`P_n = P_p = 0`.

.. _mobility:

Mobility models
^^^^^^^^^^^^^^^

The models to be used for electrons and holes can be defined in a single physical model
block or independently using two blocks. The corresponding keywords are *mobility* or
*electron_mobility* and *hole_mobility*, i.e.::

  physical_model mobility
    {
     type = ...
     ...
    }
  physical_model electron_mobility
    {
     type = ...
     ...
    }
  physical_model hole_mobility
    {
     type = ...
     ...
    }

When using the first approach, both carriers will use the same model, and param
eters provided in the input file will also be used by both carriers. When mixing the
different definitions, the blocks electron_mobility and hole_mobility will override
the common mobility block.

The default model is the constant mobility model. The parameters for the different
mobility models are needed for both electrons and holes. In the material files they are
specified with a common keyword in arrays, e.g.::

  [mobility/constant]
  # electrons holes
  mu_max = (1400.0 , 250.0)
  exponent = (1.0 , 2.1)
  
**Constant mobility model**
  
The constant mobility model (identifier **constant**) assumes a mobility which depends
only on temperature by means of the following formula:

.. math::
   :nowrap:
   :label:
   
   \begin{equation}
   \mu_{const} = \mu_0 (T/T_0)^{-\gamma}
   \end{equation}
   
| 
   
In the material data file :math:`\mu_0` and :math:`\gamma` 
 have to be specified with the keywords **mu_max** and
**exponent** . :math:`\mu_0` can be ovverridden from the **physical_model** section using the keyword
mu or from the **Region** sections using the keywords **mu_e** and **mu_h**.

**Doping dependent mobility model**

The doping dependent mobility model ( identifier **doping_dependent** ) implements two
models for mobility depending on the total doping density and the temperature. The
model that is used depends on the value of the **mobility_formula** parameter.

*Model by Masetti* et al. [Cit-D]_

The model by Masetti et al. is identified by **mobility_formula = 1** . It uses the following
formula:

.. math::
   :nowrap:
   :label:
   
   \begin{align}
   \mu = \mu_{min,1} * \mathrm{e}^{-P_c/N} +
   \frac{\mu_{const} - \mu_{min,2}}{1 + (N/C_r)^\alpha} -
   \frac{\mu_1}{1 + (C_s/N)^\beta}
   \end{align}

where N is the total doping density and :math:`\mu_const` the mobility obtained from the constant
mobility model. The parameters are specified in the material file as given in :ref:`Table64` .

*Model by Arora* [Cit-E]_

The model by Arora is identified by **mobility_formula = 2** . It reads:

.. math::
   :nowrap:
   :label:
   
   \begin{align}
   \mu = \mu_{min} +& \frac{\mu_d}{1+(N/N_0)^{A^*}} \\
   \text{with} \nonumber \\
   \mu_{min} = A_{min}(T/T_0)^{\alpha_m},& \quad \mu_d = A_d(T/T_0)^{\alpha_d} \nonumber \\
   N_0 = A_N(T/T_0)^{\alpha_N},& \quad A^* = A_a(T/T_0)^{\alpha_a} \nonumber 
   \end{align}
   
The parameters are given in :ref:`Tabl65` .

**Field dependent mobility model**

The field dependent mobility model describes the degradation of mobility at high driving
fields. It is identified by the identifier **field_dependent** . The electric field component
in directon of the current flow or the gradient of the electro-chemical potential can be
chosen as driving force:

**driving_force = efield | grad_fermi | field_parameter**

The default driving force is the gradient of the corresponding electro-chemical potential
:math:`\nabla\phi`. **field_parameter** uses a field parameter given by
:math:`\sqrt{E\cdot\nabla\phi}` as driving force [?].
The model is based on the Caughey-Thomas model, refined by Canali [Cit-F]_ :

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
   \beta = \beta_0 (T/T_0)^b
   \]

**|E|** is the modulus of the driving field, :math:`\mu_{lowfield}` is the low-field mobility. For the latter 
one can specify the model to be used using the parameter **lowfield_model** . As default
the doping dependent model is used.
There are two models for vsat, identified with **Vsat_Formula = 1 and 2** . 

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


The parameters for the field dependent mobility model are summarized in :ref:`Table66` .

.. index:: single:Boundary

Boundary conditions
^^^^^^^^^^^^^^^^^^^

Boundary conditions are implemented for ohmic contacts, Schottky contacts, free surfaces
and interfaces. Contacts are boundary models that allow a nonzero normal electrical
current. The applied voltage is specified with the option **voltage** . A variable can be
assigned to this, using the @-syntax.

For a finer control of the behaviour at electrical contacts, the options **zero_field** ,
**zero_grad_fermi_e** and **zero_grad_fermi_h** can be used, which when set to **true** will
impose zero normal electric field and zero normal gradient of the electron and hole
electro-chemical potential, respectively. The ohmic contact (identifier **ohmic**) has no
further parameters.

A Schottky contact (identifier **schottky**) has the additional parameter barrier,
which signifies the energy difference between the semiconductor band edge and the fermi
energy in the metal. As default, the barrier is taken with respect to the conduction
band. By specifying **band = v** the barrier can be imposed with respect to the valence
band (p-type contact). Alternatively, the metal work function can be defined using the
keyword **work_function** . Note, however, that its value has to be aligned with the band
energies given in the material files for the other materials. The **fixed_barrier** controls
the behaviour of the barrier height for strained materials. If it is set to true, the barrier
will be independent of strain (default behaviour). If it is set to **false** , the given barrier
is used as barrier for the unstrained case and will depend on strain during simulation. If
the metal work function is specified, the barrier will be strain dependent as default.

The type of boundary model is chosen by the parameter **type** , e.g. type = schottky.
A surface recombination can be specified for electrons and holes by using the keywords
**rec_velocity_n** and **rec_velocity_p** .

The free surface or interface model (identifier **interface**) models surface charges and
surface recombination. Two modes are possible for the surface charge model:

* :term:`constant charge` a constant charge can be assigned by specifying only the sheet carrier
density Ns in :math:`cm^{-2}`. The sheet charge density will then equal Ns multiplied by the
elementary charge e. A positive Ns produces a positive surface charge.

* :term:`electronic surface states` in this case the surface charge is produced by electrons oc-
cupying a surface state with a density of states in form of a delta function. The
density of occupied states then reads

..   math::
     :nowrap:
     :label:
     
     \[
     n_s = \frac{N_s}{1 + \frac{1}{g}\exp(\frac{E_c - \Delta E_s - e\varphi + e\phi_n}{k_BT})}
     \]

The density of states :math:`N_{s}` is specified by **Ns** , the energy of the state with respect
to the conduction band :math:`\Delta E_s` by Es. *g* denotes the multiplicity of the state and
defaults to 2. It can be changed by assigning a value to *g* .

Surface recombination is switched on by setting the keyword **surface_rec** to **true** .
The model used for surface recombination is formally the same as the bulk SRH recom-
bination. However, the only parameters are the recombination velocities, which have to
be specified (in cm/s) in the input file using the keywords **rec_vel_e** and **rec_vel_h** .

.. _PhysicSection:

Physics section
---------------

Options for controlling the drift-diffusion semiconductor models can be specified in the
Physics section. The corresponding paramaters are given in :ref:`Table67` . When **model** is
not specified, the default semiconductor model based on bulk :math:`k \cdot p` theory is used.

The **electron_quantum_density** and **hole_quantum_density** will use the particle 
densities calculated from the corresponding **quantum_density** simulation. In regions 
where no quantum density is available, the classical density will be used. The
**electron_quantum_density** and **hole_quantum_quantum_density** keywords can be used
also in the **Region** sections to be able to use different quantum density simulations in
different regions. For further options regarding selfconsistent Schrodinger-Poisson/DriftDiffusion 
calculations see :ref:`specialoptionsSchro` .

The **strain_simulation** option is used to specify the simulation that provides strain
in the case of strained systems. If it is omitted, an unstrained system is assumed for the
drift-diffusion calculation.

The **thermal_simulation** option is used to specify the simulation that provides the
lattice temperature for non-isothermal simulations. If it is omitted, the simulation tem-
perature as provided in the **Simulation** section of the input file (or, if not provided, the
default value of 300 K) is used.

Simple semiconductor model
^^^^^^^^^^^^^^^^^^^^^^^^^^

When specifying **model = simple** a very simple semiconductor model is used. For this
model one has to provide conduction and valence band edges and the effective density of
states masses in the **Region** sections. The corresponding keywords are given in :ref:`Table68` .

Default semiconductor model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default semiconductor model uses a bulk :math:`k \cdot p` model to calculate the band parameters.
It can be chosen explicitly by **model = default** . The model reads all needed parameters
from the material data file.

The band parameters are calculated considering locally strain and lattice temperature
as obtained from the corresponding simulations specified using the **strain_simulation**
and **thermal_simulation** keywords.

.. _specialoptionsSchro: 

Special options for Schrodinger-Poisson calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TIBERCAD is able to do selfconsistent Schrodinger-Poisson or Schrodinger-Drift-Diffusion
calculations. For this purpose, **electron_quantum_density** or **hole_quantum_density**
has to be specified for at least one region, and a selfconsistent simulation should be
defined in the **Selfconsistent** block. The following options { to be specified in the
**Physics** section { control the behaviour of the selfconsistent simulation.

* **use_density_predictor = bool** When set to true, a predictor-corrector scheme will
be adopted in the selfconsistent cycle. The Poisson/Drift-Diffusion solver does not
just take the particle densities as given by the Schrodinger calculation, but it will
assume a dependency of the density on the potentials of the form

..  math::
    :nowrap:
    :label:
    
    \begin{equation}
    \rho(\varphi, \phi_n, \phi_p) = \frac{\rho_{\mathrm{quantum}}(\varphi^0, \phi_n^0, \phi_p^0)}{\rho_{\mathrm{classical}}(\varphi^0, \phi_n^0, \phi_p^0)}\rho_{\mathrm{classical}}(\varphi, \phi_n, \phi_p)
    \end{equation}

where :math:`(\phi^0, \phi_n^0, \phi_p^0)` are the potentials for which the quantum density was calculated. 
**use_density_predictor = true** is the preferred method for selfconsistent
Schrodinger-Poisson/Drift-Diffusion calculations, however it is not enabled by default.

* **embracing_length = double** When the domain of the quantum simulation is smaller
than the domain of the full simulation, the boundary conditions for the Schrodinger
equation will disturb the transfer from classical to quantum density. By defining an
embracing region of a certain extension (specified in meters), a gradual transition
from classical to quantum density will be done instead of an abrupt one, using as
effective density :math:`\rho = x \cdot \rho_{quantum} + (1 - x) \cdot \rho_{classical}`. The default is no embracing
region at all (zero extension).

* **cutoff = double** if an embracing region is used, a part of this region near the boundary
of the quantum region can be cut off so that only the classical density is considered
in that part. **cutoff** is specified as a percentage of the embracing length and should
therefore be between 0.0 and 1.0.

* **plot_embracing_regions = bool** Whereas the automatic creation of the embracing
region in 1D is a very simple task, it is a more difficult one in higher dimensions.
By setting this flag to **true** , 
the embracing region and the mixing coefficient x will
be plotted for a visual control of the quality of the embracing region. The default
is **false** .

Solver section
--------------

Many of parameters for the numerical solver depend on the type of solver being used
and on the device to be simulated. :ref:`Table69` lists the options that are independent on
the type of solver used.
The linear and nonlinear solvers to be used can be chosen using the keywords **linear_solver**
and **nonlinear_solver** , respectively. For the nonlinear solver one can chose between
the PETSc implementation ( **petsc** ) and the TIBERCAD implementation ( **tiber** ) of line
search. When using the TIBERCAD nonlinear solver, one can additionally chose between
the PETSc ( **petsc** ) or PARDISO ( **pardiso** ) linear solvers. The possible combinations are::

  nonlinear_solver = petsc
  or
  nonlinear_solver = tiber
  linear_solver = petsc | pardiso

Parameters for PETSc solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`Table610` and :ref:`Table611` list all solver parameters significant for the PETSc linear and
nonlinear solvers. A more detailed description of the most important parameters follows.
The **ksp_type** specifies the type of Krylov subspace method to be used. The mostly
used methods are:

* **bcgs** A stabilized version of the biconjugate gradient method. This one works better in
1D than **bcgsl**.

* **bcgsl** (default) A modified version of **bcgs**.

* **gmres** Generalized minimal residual method.

The **pc_type** specifies the type of preconditioner to be used. The most useful ones
are:

* **ilu** (default) Incomplete LU factorization. Does not work for materials with high band
gap.

* **jacobi** Jacobi preconditioning (diagonal scaling).

* **composite** Combination of **ilu** and **jacobi** .

The **ls_max_step** parameter defines an upper bound of the l2-norm of the nonlinear
line search step. It should be not too big to prevent the algorithm from diverging, but
also not too small to minimize the number of iterations. Values between 1 and 10 should
be a good choice.

The **nonlin_step_tol** defines at which line search step size (in :math:`l_2-norm`) the algorithm
stops, i.e. assumes to have reached convergence. **nonlin_step_tol** is measured in eV.

Parameters for the TIBERCAD nonlinear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`Table12` summarizes the parameters used for the TIBERCAD implementation of the
line search algorithm.

The stopping criterion based on the line search step uses the maximum norm of the
nonlinear step, i.e. convergence is controlled locally. In addition to the parameters in
:ref:`Table6-12` one has to provide also parameters for the linear solver.

Parameters for the PARDISO linear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE:: |warn| 
    the PARDISO linear solver is currently not included in the distribution.
    Currently the PARDISO interface has no adjustable parameters.
    
----

.. _Table61: 

   Table 6.1

+---------------------+----------------------------------------------------+------------+
|Key                  |Description                                         | Units      |
+---------------------+----------------------------------------------------+------------+
|Ec                   | Conduction band edge                               |  eV        |
+---------------------+----------------------------------------------------+------------+
|Ev                   | Valence band edge                                  |  eV        |
+---------------------+----------------------------------------------------+------------+
|eQFermi              | Electro-chemical potential of electrons            |  eV |61.A| |
+---------------------+----------------------------------------------------+------------+
|hQFermi              | Electro-chemical potential of holes                |  eV |61.B| |
+---------------------+----------------------------------------------------+------------+
|Ec0                  | Conduction band edge without electric potential    |  eV        |
+---------------------+----------------------------------------------------+------------+
|Ev0                  | Valence band edge without electric potential       |  eV        |
+---------------------+----------------------------------------------------+------------+
|Eg                   | Band gap                                           |  eV        |
+---------------------+----------------------------------------------------+------------+
|ConductionBands      | Minima of all conduction bands                     |  eV        |
+---------------------+----------------------------------------------------+------------+
|ValenceBands         | Maxima of all valence bands                        |  eV        |
+---------------------+----------------------------------------------------+------------+
|ElPotential          | Electric potential                                 |  V         |
+---------------------+----------------------------------------------------+------------+
|eDensity             | Electron density                                   |  cm-3      |
+---------------------+----------------------------------------------------+------------+
|hDensity             | Hole density                                       |  cm-3      |
+---------------------+----------------------------------------------------+------------+
|eMobility            | Electron mobility                                  |  cm2V-1s-1 |
+---------------------+----------------------------------------------------+------------+
|hMobility            | Hole mobility                                      |  cm2V-1s-1 |
+---------------------+----------------------------------------------------+------------+
|eConductivity        | Electron conductivity                              |  S/cm      |
+---------------------+----------------------------------------------------+------------+
|hConductivity        | Hole conductivity                                  |  S/cm      |
+---------------------+----------------------------------------------------+------------+
|ElField              | Electric Field                                     |  Vcm-1     |
+---------------------+----------------------------------------------------+------------+
|Polarization         | Electric Polarization                              |  Cm-2      |
+---------------------+----------------------------------------------------+------------+
|CurrentDensity       | Total current density                              |  Acm-2     |
+---------------------+----------------------------------------------------+------------+
|eCurrentDensity      | Electron current density                           |  Acm-2     |
+---------------------+----------------------------------------------------+------------+
|hCurrentDensity      | Hole current density                               |  Acm-2     |
+---------------------+----------------------------------------------------+------------+
|IonizedDonors        | Ionized donor density                              |  cm-3      |
+---------------------+----------------------------------------------------+------------+
|IonizedAcceptors     | Ionized acceptor density                           |  cm-3      |
+---------------------+----------------------------------------------------+------------+
|IonizedElectronTraps | Ionized electron trap density                      |  cm-3      |
+---------------------+----------------------------------------------------+------------+
|IonizedHoleTraps     | Ionized hole trap density                          |  cm-3      |
+---------------------+----------------------------------------------------+------------+




+---------------------+----------------------------------------------------+------------+
|Key                  |Description                                         | Units      |
+---------------------+----------------------------------------------------+------------+
|eThElPower           | Electron thermoelectric power                      |  V K-1     |
+---------------------+----------------------------------------------------+------------+
|hThElPower           | Hole thermoelectric power                          |  V K-1     |
+---------------------+----------------------------------------------------+------------+
|NetRecombination     | The net recombination rate for each recombination  |  cm-3s-1   |
+---------------------+----------------------------------------------------+------------+
|                     |  generation model and the total rate               |            |
+---------------------+----------------------------------------------------+------------+
|ContactCurrents      | The electric current on each contact               |  A/cm3-d   |
+---------------------+----------------------------------------------------+------------+

Table 6.1: Solution/Plot variables

| 
| 


.. _Table62: 

   Table 6.2

+--------------------+---------+
|parameter           | keyword |
+--------------------+---------+
|:math:`\tau_min`    |  taumin | 
+--------------------+---------+
|:math:`\tau_max`    |  taumax | 
+--------------------+---------+
|:math:`N_ref`       |  Nref   | 
+--------------------+---------+
|:math:`\gamma`      |  gamma  | 
+--------------------+---------+
|:math:`E*`          |  Etrap  | 
+--------------------+---------+
|:math:`\alpha`      |  Talpha | 
+--------------------+---------+
|:math:`\beta`       |  Tcoeff | 
+--------------------+---------+

Table 6.2: SRH material data file parameters

| 
| 


.. _Table63: 

   Table 6.3

+--------------------+---------+
|parameter           | keyword |
+--------------------+---------+
|:math:`\tau_n`      |   tau_n |
+--------------------+---------+
|:math:`\tau_p`      |   tau_p |
+--------------------+---------+
|:math:`E*`          |   E_t   |
+--------------------+---------+

Table 6.3: SRH input file parameters

| 
| 


.. _Table64: 

   Table 6.4

+--------------------+---------+
|parameter           | keyword |
+--------------------+---------+
|:math:`\mu_{min,1}` |  mumin1 |
+--------------------+---------+
|:math:`\mu_{min,2}` |  mumin2 |
+--------------------+---------+
|:math:`\mu_1`       |  mu1    |
+--------------------+---------+
|:math:`P_c`         |  Pc     |
+--------------------+---------+
|:math:`C_r`         |  Cr     |
+--------------------+---------+
|:math:`C_s`         |  Cs     |
+--------------------+---------+
|:math:`\alpha`      |  alpha  |
+--------------------+---------+
|:math:`\beta`       |  beta   |
+--------------------+---------+

Table 6.4: Data file parameters for the mobility model by Masetti et al.

| 
| 


.. _Table65: 

   Table 6.5

+-------------------------+---------+
|parameter                | keyword |
+-------------------------+---------+
|:math:`A_min`            |  mumin  |
+-------------------------+---------+
|:math:`A_d`              |  mud    |
+-------------------------+---------+
|:math:`A_N`              |  N0     |
+-------------------------+---------+
|:math:`A_a`              |  A      |
+-------------------------+---------+
|:math:`\alpha_m`         |  am     |
+-------------------------+---------+
|:math:`\alpha_d`         |  ad     |
+-------------------------+---------+
|:math:`\alpha_N`         |  aN     |
+-------------------------+---------+
|:math:`\alpha_a`         |  aA     |
+-------------------------+---------+

Table 6.5: Data file parameters for the mobility model by Arora.

| 
| 


.. _Table66: 

   Table 6.6

+-------------------------+---------+
| parameter               | keyword |
+-------------------------+---------+
|:math:`\beta_0`          | beta0   |
+-------------------------+---------+
|:math:`b`                | betaexp |
+-------------------------+---------+
|:math:`v_{sat,0}`        | vsat0   |
+-------------------------+---------+
|:math:`\gamma`           | vsatexp |
+-------------------------+---------+
|:math:`A_vsat`           | A_vsat  |
+-------------------------+---------+
|:math:`B_vsat`           | B_vsat  |
+-------------------------+---------+
|:math:`v_min`            | vsat_min|
+-------------------------+---------+

Table 6.6: Data file parameters for the mobility model by Arora.

| 
| 


.. _Table67: 

   Table 6.7

+--------------------------+-------------------------------+------------------------------------------------------------+
| keyword                  | possible values               | description                                                |
+--------------------------+-------------------------------+------------------------------------------------------------+
| model                    | (see following subsections)   | the model to use for the description of the conduction and |
+--------------------------+-------------------------------+------------------------------------------------------------+
|                          |                               |      valence band properties                               |
+--------------------------+-------------------------------+------------------------------------------------------------+
| statistics               | B | FD                        | Boltzmann (default) or FermiDirac statistics               |
+--------------------------+-------------------------------+------------------------------------------------------------+
| strain_simulation        | name                          | the strain simulation to be used                           |
+--------------------------+-------------------------------+------------------------------------------------------------+
| thermal_simulation       | name                          | the thermal simulation to be used                          |
+--------------------------+-------------------------------+------------------------------------------------------------+
| electron_quantum_density | name                          | the quantum density simulation to be used for the electron |
+--------------------------+-------------------------------+------------------------------------------------------------+
|                          |                               |     density                                                |
+--------------------------+-------------------------------+------------------------------------------------------------+
| hole_quantum_density     | name                          | the quantum density simulation to be used for the hole     |
+--------------------------+-------------------------------+------------------------------------------------------------+
|                          |                               |     density                                                | 
+--------------------------+-------------------------------+------------------------------------------------------------+

Table 6.7: Common options for the drift-diffusion semiconductor models

| 
| 


.. _Table68: 

   Table 6.8

+--------------------+-----------------------------------------+
| keyword            | description                             |
+--------------------+-----------------------------------------+
| Ec                 | conduction band edge (eV)               |
+--------------------+-----------------------------------------+
| Ev                 | valence band edge (eV)                  |
+--------------------+-----------------------------------------+
| m_dos_e            | conduction band effective DOS mass (me) |
+--------------------+-----------------------------------------+
| m_dos_h            | valence band effective DOS mass (me)    |
+--------------------+-----------------------------------------+

Table 6.8: Parameters for the simple semiconductor model

| 
| 


.. _Table69: 

   Table 6.9

+----------------------------+------------------------------------------------------------------------------------------+
| keyword                    | description                                                                              |
+----------------------------+------------------------------------------------------------------------------------------+
| coupling                   | defines which equations to couple together. poisson: solve only poisson eq., electrons:  |
+----------------------------+------------------------------------------------------------------------------------------+
|                            |    electrons and poisson,                                                                |
+----------------------------+------------------------------------------------------------------------------------------+
|                            |   holes: holes and poisson, current: only electron and hole currents, <full>: the fully  |
+----------------------------+------------------------------------------------------------------------------------------+
|                            |    coupled system                                                                        |
+----------------------------+------------------------------------------------------------------------------------------+
| el_qfermi_level            | the spatially constant electrochemical potential for the electrons (for quasi-equilibrium|
+----------------------------+------------------------------------------------------------------------------------------+
|                            |     calculations)                                                                        |
+----------------------------+------------------------------------------------------------------------------------------+
| hl_qfermi_level            | the spatially constant electrochemical potential for the holes (for quasi-equilibrium    |
+----------------------------+------------------------------------------------------------------------------------------+
|                            |    calculations)                                                                         |
+----------------------------+------------------------------------------------------------------------------------------+
| integration_order          | order of the numerical gauss integration. Default is 2                                   |
+----------------------------+------------------------------------------------------------------------------------------+
| current_integration_method | method for the calculation of the contact cutrrents. surfint: integrate the local current| 
+----------------------------+------------------------------------------------------------------------------------------+
|                            |    density over                                                                          |
+----------------------------+------------------------------------------------------------------------------------------+
|                            |   the contact surface, <rstf>: use the Ramo-Shockley test functions, gives better results|
+----------------------------+------------------------------------------------------------------------------------------+
| local_scaling              | apply a local scaling scheme which leads to better scaled matrices. <true> or false      |
+----------------------------+------------------------------------------------------------------------------------------+
| exact_newton               | use exact or approximate (without some parts in the jacobian) Newton. <true> or false    |
+----------------------------+------------------------------------------------------------------------------------------+


Table 6.9: Solver independent parameters [#]_

| 
| 

.. _Table610: 

   Table 6.10

+------------------+--------------------------------------------------------------------------------------+-------------+
| keyword          | description                                                                          |  default    |
+------------------+--------------------------------------------------------------------------------------+-------------+
| nonlin_rel_tol   |    relative tolerance for the residual l2-norm(with respect to first nonlinear step) |  10e-9      |
+------------------+--------------------------------------------------------------------------------------+-------------+
| nonlin_abs_tol   |    absolute tolerance for the residual l2-norm                                       |  10e-15     |
+------------------+--------------------------------------------------------------------------------------+-------------+
| nonlin_step_tol  |    tolerance for the l2-norm of the nonlinear step                                   |  10e-3      |
+------------------+--------------------------------------------------------------------------------------+-------------+
| ls_max_step      |    the maximum l2-norm for the line searchstep, in units of eV                       |  1          |
+------------------+--------------------------------------------------------------------------------------+-------------+
| nonlin_max_it    |    maximum number of nonlinear iterations                                            |  20         |
+------------------+--------------------------------------------------------------------------------------+-------------+

Table 6.10: Parameters for the PETSc nonlinear solver

| 
| 


.. _Table611: 

   Table 6.11

+-----------------------------+---------------------------------------------+-------------+
| keyword                     | description                                 | default     |
+-----------------------------+---------------------------------------------+-------------+
| ksp_type                    | the linear solver type                      | bcgsl       |
+-----------------------------+---------------------------------------------+-------------+
| pc_type                     | the preconditioner type                     | ilu         |
+-----------------------------+---------------------------------------------+-------------+
| lin_rel_tol                 | relative tolerance for the linear solver    | 1e-6        |
+-----------------------------+---------------------------------------------+-------------+
| lin_abs_tol                 | absolute tolerance for the linear solver    | 10e-50      |
+-----------------------------+---------------------------------------------+-------------+
| lin_max_it                  | maximum number of linear iterations         | 500         |
+-----------------------------+---------------------------------------------+-------------+


Table 6.11: Parameters for the PETSc linear solver [#]_

| 
| 


.. _Table612: 

   Table 6.12

+------------------+-----------------------------------------------------------------------------------+---------+
| keyword          | description                                                                       | default |
+------------------+-----------------------------------------------------------------------------------+---------+
| nonlin_rel_tol   | relative tolerance for the residual l2-norm(with respect to first nonlinear step) | 10e-9   |
+------------------+-----------------------------------------------------------------------------------+---------+
| nonlin_abs_tol   | absolute tolerance for the residual l2-norm                                       | 10e-15  |
+------------------+-----------------------------------------------------------------------------------+---------+
| nonlin_step_tol  | tolerance for the maximum norm of the nonlinear step (eV)                         | 10e-3   |
+------------------+-----------------------------------------------------------------------------------+---------+
| nonlin_max_it    | maximum number of nonlinear iterations                                            |  20     |
+------------------+-----------------------------------------------------------------------------------+---------+

Table 6.12: Parameters for the TIBERCAD line search




.. |more| image:: more.*
    :scale: 50%

.. |warn| image:: warn.*
    :scale: 50%

.. |idea| image:: idea.*
    :scale: 50%

.. |61.A| replace:: :math:`\phi_n`

.. |61.B| replace:: :math:`\phi_p`

.. _fastlink: http://www.tibercad.org/documentation/tutorial/

.. rubric:: Footnotes

.. [#] the Default value is given in brackets.
.. [#] the linear tolerance gets automatically decreased after each nonlinear step.
