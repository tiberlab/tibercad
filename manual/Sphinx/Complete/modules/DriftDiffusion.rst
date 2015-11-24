
.. _DriftDiffusionTheory:


Drift-diffusion simulation of electrons and holes
=================================================

Theory
----------------


The semi-classical transport simulation of electrons and holes is based on the drift-diffusion
approximation (see [Selberherr]_).


Beside the electric potential the electro-chemical potentials are used as variables such
that the system of PDEs to be solved reads as follows

.. math::
   :label: dd_eq_ddsystem
   
   -\nabla(\varepsilon\nabla\varphi - \mathbf{P}) & =  -e(n - p - N_d^+ + N_a^-) \\
   -\nabla(\mu_n n ( \nabla\phi_n + P_n \nabla T)  ) & =  R \\
   -\nabla(\mu_p p (\nabla\phi_p + P_p \nabla T) ) & =  -R 
 
:math:`P` is the electric polarization due to e.g. piezoelectric effects.
:math:`N_d^+` and :math:`N_a^-` are the densities of ionized donors and acceptors, respectively.
:math:`R` is the net recombination rate, i.e. recombination rate minus generation rate, and :math:`P_n` and :math:`P_p` are the electron
and hole thermoelectric powers, respectively. The models for the mobilities and the net
recombination rates can be specified in the ``Physics`` section as described in the
following.

..
  ..  index:: double:DriftDiffusion;Solution

Solution/Plot variables
-----------------------

The solution variables available for plotting and for interaction with other modules are
given in :ref:`Plotting variables<dd_solutions>` .


Module options
---------------------

The following options influence the behaviour of the Drift-Diffusion module:

 ``coupling`` : string
    defines which equations to solve. 
    The default is ``full``, meaning that
    the full system consisting of the Poisson, electron continuity and hole continuity
    equations is solved. Other possible values are ``poisson`` (for equilibrium calculations),
    ``electrons`` or ``holes``. For the last two cases local equilibrium is assumed such that
    :math:`\phi_n  =  \phi_p` .
 
 ``enforce_local_charge_neutrality`` : bool 
    If set to true, local charge neutrality will
    be enforced by setting the charge density to zero. This may be useful for solving
    the Poisson equation involving only dielectrics.

 ``guess_el_qfermi`` : double 
     If this option is set, then before resolving the system the
     given number will be set as a guess for the electron electro-chemical potential.

 ``guess_hl_qfermi`` : double 
     If this option is set, then before resolving the system the
     given number will be set as a guess for the hole electro-chemical potential.

 ``default_boundary_condition`` : string 
        With this option the user can control the default 
        boundary condition for the electric field on all external boundaries without explicit boundary model. 
        Possible values are zero field (default), or zero displacement.
        The two differ only in presence of electric polarization fields.

 ``quadrature_rule`` : string 
       This option allows to chose between trapezoidal and Gauss
       type numeric integration rules. The default rule is gauss, but in some cases trapez
       may prevent density peaks near badly resolved material interfaces.

 ``save_state`` : bool
        If set to ``true`` the current solution will be written to a compressed 
        file after each solve. The file name follows the same rules as the result files,
        having file extension ``.tsv`` .

 ``load_state`` : string 
       Reload a formerly saved solution. The provided string value has to be the relative or
       absolute path to a ``.tsv`` file created using the option ``save_state``.

 ``solve_after_load`` : bool
       If set to ``true``, the system will be solved after having
       reloaded a saved state. Otherwise it will not be solved, which is the default behaviour.

.. warning::  
            Currently the reload of saved solutions *only works correctly when using the identical mesh*. 
            Otherwise there will be undefined behaviour or failure. 


Solver section
--------------------

The ``Solver`` section of the Drift-Diffusion module refers to a nonlinear solver. 
See  :ref:`Nonlinear_solver` for details on nonlinear solver options.

Physics section
--------------------

The ``Physics`` block contains generic options for the bulk physical model and the definition
of submodels. The generic options are:

 ``thermal_simulation`` 
     If you are doing coupled electrothermal simulations,
     you have to specify the name of the thermal simulation providing the lattice temperature.

 ``strain_simulation`` 
     If you are doing simulations on strained systems, you have to
     specify the name of a strain simulation. The strain values obtained from
     this model will be used to calculate strain dependent parameters like
     band parameters and piezoelectric polarizations.

 ``relax_polarization`` 
     With this option one can specify a global relaxation factor
     for the electric polarization field. This can be useful if the amount of total electric
     polarization has to be treated as fitting parameter.

In the following we describe all submodels. As mentioned in the Introduction (REF HERE)
submodels can be restricted to a subset of simulation regions.


Band parameter models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. The model and parameters to calculate band parameters like band edge, effective density of states and particle densities

The model and parameters to calculate band parameters like band edges, band gaps etc. 
can be controlled by special submodel blocks in different ways:

1. a single ``band_properties`` block containing parameters for conduction and for valence band

2. ``conduction_band`` and ``valence_band`` blocks to control each band independently

Band parameter models may be defined through the block *density_of_states*, described in the following.


Density of states (DOS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The density of states for electrons and/or holes can be set using the ``density_of_states`` submodel block inside the ``band_properties`` block, as follows::

   band_properties # or conduction_band or valence_band
   {
     density_of_states [type]
     {
       ...
     }
   }

There are several types of density of states models.


3D Bulk DOS (default)
............................................

This is the default model describing the 3D DOS of a parabolic band semiconductor.
Band edge energies and the DOS mass are taken from the database, or can be provided from the input file.
In the current version, if band  parameters have to be given from the input file, it is necessary to use the ``conduction_band`` and ``valence_band`` blocks instead of ``band_properties``.
In the ``density_of_states`` subblock, one can then give the band edge and the DOS mass the keywords ``level``, and ``dos_mass``, respectively.
A degeneracy can be specified with the ``degeneracy`` keyword.

Bulk kp DOS
............................................
This model for bulk band parameters can be chosen by specifying ``density_of_states bulk_kp``. It calculates the band edge energies and masses from bulk kp theory, including Pikus-Bir strain corrections.
To include strain corrections, the keyword ``strain_simulation`` has to be used, providing the name of the module instance which calculates strain.
In the current version it is not possible to provide kp parameters from the input file.
For example ::

  band_properties 
    {
      density_of_states bulk_kp 
      {
        strain_simulation = strain
      }
    }


Quantum DOS
............................................

With this model, the local DOS is obtained from the solution of the :math:`Schr\ddot{o}dinger` equation in a defined system. The following keywords are used to control this model :


 ``quantum_simulation`` 
     The name of a quantum density simulation. This will use the quantum mechanical particle density in the regions it was calculated. More than one simulations can be specified as a vector. In this case the sum of all quantum densities is taken.
If  ``quantum_simulation`` is specified, the following additional options control the mixing between
classical and quantum density:

  ``barrier_regions``
     Regions which can be regarded as pure barriers. In these regions a classical density will be added using the barrier  materials bulk band edge. 

A further submodel block may be added to define an alternative  model to be used to obtain the DOS inside the barrier regions defined by *barrier_regions*, or when the quantum density is not available, for example because it has not been yet calculated in the simulation. 
This block is called *classical_DOS* , and its type may be any one of the types defined for *density_of_states*  model, as follows::

   classical_DOS [type]
     {
       ...
     }
   


For example: ::

  conduction_band
    {
      density_of_states quantum
      {
        # where to get the quantum density from
        quantum_simulation = quantum_el
        barrier_regions = buffer_quantum

        classical_DOS bulk_kp
        {
          strain_simulation = strain
        }

      }
    }

Here, quantum DOS, and thus quantum density, for electrons, is obtained from the simulation *quantum_el*, for example an EFA calculation. A classical DOS is defined, to be obtained with the model *bulk_kp*.


Constant DOS
............................................

The constant DOS model is identified by the identifier ``constant``. The constant DOS has the following form, for :math:`E_{min} \le E \le E_{max}`:

.. math::
   g(E) dE = \frac{N_0}{E_{max} - E_{min}} dE

and :math:`g(E) = 0` otherwise. In the input file the constant DOS is set using three keywords, ``N0``, ``level`` and ``Ewidth``. The first one is used to set the maximum occupied density :math:`N_0` (:math:`cm^{-3}`). If a value for ``level`` is not specified in the input file, it will be set by default as the conduction band edge (:math:`E_C`) for electrons, and as the valence band edge (:math:`E_V`) for holes. The values of :math:`E_C` and :math:`E_V` are taken from material files. For electrons :math:`E_{min} =` ``level``, and :math:`E_{max} =` ``level`` + ``Ewidth``; for holes :math:`E_{max} =` ``level``, and :math:`E_{min} =` ``level`` - ``Ewidth``.



Gaussian DOS
............................................

The gaussian DOS model is identified by the identifier ``gauss``. The gaussian DOS has the following form

.. math::
   g(E, \sigma) dE = \frac{N_0}{\sigma \sqrt{2 \pi}} \exp \left( - \frac{1}{2} \frac{\left(E- E_{C,V}\right)^2}{\sigma^2} \right) dE

where :math:`N_0` (:math:`cm^{-3}`) is the maximum occupied density of carriers; :math:`\sigma` (:math:`eV`) is the variance of the gaussian, hence it is a measure of the energy disorder; :math:`E_C` and :math:`E_V` are the conduction and valence band edges respectively. Band edges are set in the material files of TiberCAD, but a different center for the gaussian can be specified with the keyword ``level = [value]``. :math:`N_0` and :math:`\sigma` can be set with the keywords ``N0 = [value]`` and ``sigma = value`` respectively.

The gaussian DOS can be used also to define a distribution of trap states (see section :ref:`DD_trapmodels`), for example::

   trap eNeutral
   {
     regions = ...
     Nt = 1e18
     Et = 1.0
     reference = cb

     density_of_states
     {
       type = gauss
       N0 = 1
       sigma = 0.1
     }
   }

According with the conventions used in TiberCAD the maximum density of traps is obtained as the product ``Nt`` :math:`\times` ``N0``, hence one can set indifferently ``Nt`` or ``N0``, leaving one parameter set to 1. The center of the gaussian is set with the value ``Et``, thus it is not necessary to provide a value for ``level``. ``Et`` is specified with respect to the conduction band (``reference = cb``), or the valence band (``reference = vb``) or the midgap (``reference = m``). The example sets a gaussian distribution of ``eNeutral`` traps whose center is 1 :math:`eV` below the conduction band.



The gaussian DOS can be used also to define a distribution of trap states (see section :ref:`DD_trapmodels`). This is particularly useful when modeling metal/semiconductor interfaces using the IDIS model (Induced Density of Interface States, see [Santoni]_). According to the IDIS model a density :math:`N_{trap}` of traps is distributed on the interface and a Charge Neutrality Level (CNL) is defined in such a way that when the fermi level is below the CNL the interface is positively charged, and when the fermi level is above the CNL the interface is negatively charged. This situation can be simulated using ``eNeutral`` traps over the CNL, and ``hNeutral`` traps below the CNL. For :math:`Alq_3` the CNL is 1.68 :math:`eV` above the valence band and :math:`N_{trap} = 2.63\times10^{19} \ cm^{-3} eV^{-1}`. Thus if we want to simulate a constant distribution of ``eNeutral`` traps above the CNL extending up to the conduction band edge (bandgap for :math:`Alq_3` is 2.9 :math:`eV`) we have to set::

  trap eNeutral
  {
    regions = ...
    Nt = 3.21e19
    Et = 1.68  #level of the constant DOS is automatically 
                fixed by this value
    reference = vb #Et is set w. r. t. valence band edge

    density_of_states
    {
      type = constant
      Ewidth = 1.22 #Alq3 bandgap - Et
      N0 = 1
    }
  }

And to set a constant distribution of ``hNeutral`` traps between the valence band edge and the CNL we write::

  trap hNeutral
  {
    regions = ...
    Nt = 4.42e19
    Et = 1.68  
    reference = vb 

    density_of_states
    {
      type = constant
      Ewidth = 1.68 #it is the same value of Et because
                     we have set Et as the CNL w. r. t. valence band
      N0 = 1
    }
  }

According with the conventions used in TiberCAD the maximum density of traps is obtained as the product ``Nt`` :math:`\times` ``N0``, hence one can set indifferently ``Nt`` or ``N0``, leaving one parameter set to 1. We have also to keep in mind that the constant distribution is normalized in energy, so ``Nt`` (or ``N0``) is automatically divided by ``Ewidth``. But the density of states :math:`N_{trap}` is already in :math:`cm^{-3} eV^{-1}` units, so we have to multiply it by 2.9 - 1.68 = 1.22 :math:`eV` for ``eNeutral`` and by 1.68 :math:`eV` for ``hNeutral`` in order to get the correct volume density for ``Nt`` (or ``N0``).




..  ``add_continuum_in_well``
     If this option is set to true, the energy level of the first state after the ones considered for the quantum
     density is used as an effective bulk band edge and a classical carrier density will be added accordingly.

.. If a quantum density is used, then it is useful to define also an embracing region
.. where the model gradually switches from a fully classical to a fully quantum density.
.. The options for the embracing are specified in a block with keyword ``embracing`` . It
.. accepts the following options:

..   ``embracing_length = double`` 
       When the domain of the quantum simulation is smaller
       than the domain of the full simulation, the boundary conditions for the :math:`Schr\ddot{o}dinger` 
       equation will disturb the transfer from classical to quantum density. By defining an
       embracing region of a certain extension (specified in meters), a gradual transition
       from classical to quantum density will be done instead of an abrupt one, using as
       effective density :math:`\rho = x\cdot\rho_{\mathrm{quantum}} + (1-x)\cdot\rho_{\mathrm{classical}}` . 
       The default is no embracing region at all (zero extension).

..   ``cutoff = double`` 
       If an embracing region is used, a part of this region near the boundary
       of the quantum region can be cut off so that only the classical density is considered
       in that part. ``cutoff`` is specified as a percentage of the embracing length and should
       therefore be between 0.0 and 1.0.

..   ``plot_embracing_region = bool`` 
       Whereas the automatic creation of the embracing region 
       in 1D is a very simple task, it is a more difficult one in higher dimensions. 
       By setting this flag to true, the embracing region and the mixing coefficient x will be
       plotted for a visual control of the quality of the embracing region. 
       The default is ``false`` .


.. _DD_recombinationmodels:

Recombination/generation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the currently available generation/recombination models. Note
that all recombination models can be applied also to surfaces/interfaces.
Recombination models are controlled by means of ``recombination`` submodel blocks
inside the Physics section. Different recombination models having the same (or no)
options can be enabled in a single statement by writing::

  recombination (model1, model2, ...) {}



Shockley-Read-Hall (SRH) recombination
............................................


The SRH recombination model can be enabled by defining a ``recombination`` submodel
of type ``srh`` .

SRH recombination is defined as follows:

.. math::
    :label: dd_eq_recsrh
    
    R_{SRH} = \frac{np - n_i^2}{(n + n_i e^{E^*/k_BT})\tau_p +
    (p + n_i e^{-E^*/k_BT})\tau_n}


:math:`E^* = E_{trap} - (E_c + E_v)/2` is the trap level with respect to the midband energy. 

:math:`n_i` is the intrinsic carrier density, :math:`\tau_n` and :math:`\tau_p` are the recombination times. 

The parameters are taken from the material database. The recombination times are dependent on temperature and
doping density, e.g.

.. math::
    
    \tau_n & =  \tau_n^0 \left(\frac{T}{T_0}\right)^{\alpha_n} e^{\beta(T/T0 - 1)} \\
    \tau_n^0 & =  \tau_{min,n} + \frac{\tau_{max,n} - \tau_{min,n}}{1 + (N/N_{ref})^\gamma}

where :math:`T_0` is the reference temperature (300 K). Table 2.3 shows the corresponding parameters 
for the material data files. The parameters for holes and electrons have to be
specified in an array, e.g. :math:`\tau_{min} = (1e-5, 3e-6)`

    

The recombination times and trap level can be overridden from the input file by using
the keywords of Table 2.4.



The SRH recombination model can be applied also to surfaces and interfaces. In this
case, you can provide the recombination velocities using the keywords ``rec_velocity_n``
and ``rec_velocity_p`` instead of ``tau_n`` and ``tau_p`` .

Direct (radiative) recombination
.................................


The direct recombination model can be enabled in the input file by by defining a
``recombination`` submodel of type ``direct`` .

Direct recombination is modeled as follows:

.. math::
   :label: dd_eq_recdirect

    R_{direct} = C(np - n_i^2)

The material data file and the input file use the same keyword C for the parameter C. The
database value can be overridden from the input file as described for SRH recombination.

Langevin (radiative) recombination
.................................

The Langevin recombination model can be enabled in the input file by defining a
``recombination`` submodel of type ``langevin`` .

Langevin recombination is modeled as follows:

.. math::
   :label: dd_eq_reclangevin

    R_{Langevin} = \gamma \frac{e}{\varepsilon} (\mu_n + \mu_p) (np - n_0 p_0)

:math:`\mu_n` and :math:`\mu_p` are electron and hole mobilities respectively and their values are taken from the proper model specified in the simulation (see :ref:`Mobility_models`); :math:`n_0` and :math:`p_0` are equilibrium electron and hole densities respectively; :math:`\gamma` is a constant that can be set using the keyword ``gamma = [value]`` (by default ``value`` is set to 1).


Auger recombination
.................................


The Auger recombination model can be enabled in the input file by defining a recombination
submodel of type ``auger`` .

Auger recombination is modeled by the following equation

.. math::
   :label: dd_eq_recauger

    R_{auger} = (C_nn + C_pp)(np - n_i^2)

with temperature dependent parameters

.. math::

    C_{\{n,p\}} = \left(A + B\frac{T}{T_0} + C\left(\frac{T}{T_0}\right)^2\right)\left(1 + H e^{-\{n,p\}/N_0}\right)
 
The parameters A;B;C;H and :math:`N_0` are taken exclusively from the database. They are
different for :math:`C_n` and :math:`C_p` and have to be specified as arrays with keywords ``A``, ``B``, ``C``, ``H``, ``N0``,
e.g. ``A = (1e-31, 1e-32)``. The calculated values for :math:`C_n` and :math:`C_p` can be overridden from
the input file by specifying values for the keywords ``Cn`` and ``Cp`` .


Optical generation
.................................


A very simple model for photoelectric generation of electron-hole pairs is implemented
in tiberCAD. It is enabled by specifying a ``generation`` submodel of type optical The
model imposes a constant generation rate which has to be provided by the keyword G in
units of :math:`(\mathrm{cm}\cdot\mathrm{s})^{-1}` . 

.. note:: 
          Usually the simulation should define a sweep on the value
          of G from 0 to the desired generation.


.. _ThermoelectricPower:

Thermoelectric power models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermoelectric power models are the same for electrons and holes. 

The keyword is  ``thermoelectric_power`` , i.e ::

  thermoelectric_power [type]
  {
  }
    
The type can be ``constant`` (i.e. the thermoelectric powers are read from the
database) or ``diffusivity_model`` where the thermoelectric powers are computed by

.. math::
   :label: dd_eq_thermopower

    P_n & = - \frac{k_b}{q}\left( \frac{5}{2} + \frac{e \phi_n + E_c - e \varphi}{k_b T} \right) \\
    P_p & = \frac{k_b}{q}\left( \frac{5}{2} - \frac{e \phi_p + E_v - e \varphi}{k_b T} \right)

The default is :math:`P_n = P_p = 0`

.. _Mobility_models:

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
............................................




The constant mobility model (identifier ``constant`` ) assumes a mobility which depends
only on temperature by means of the following formula:

.. math::
   :label: dd_eq_muconst

    \mu_{const} = \mu_0 (T/T_0)^{-\gamma}

In the material data file :math:`\mu_0` and :math:`\gamma`  
have to be specified with the keywords ``mu_max`` and ``exponent``. 
:math:`\mu_0` can be overridden from the ``physical_model`` section using the keyword
``mu`` or from the Region sections using the keywords ``mu_e`` and ``mu_h`` .


Doping dependent mobility model
............................................


The doping dependent mobility model (identifier ``doping_dependent`` ) implements two
models for mobility depending on the total doping density and the temperature. The
model that is used depends on the value of the ``mobility_formula`` parameter.

**Model by Masetti et al.**

The model by Masetti et al.  Masetti is identified by ``mobility_formula`` = 1. It uses the following
formula:

.. math::
   :label: dd_eq_dopdep1
    
    \mu = \mu_{min,1} \cdot \mathrm{e}^{-P_c/N} + \frac{\mu_{const} - \mu_{min,2}}{1 + (N/C_r)^\alpha} - \frac{\mu_1}{1 + (C_s/N)^\beta}


where N is the total doping density and :math:`\mu_{const}` the mobility obtained from the constant
mobility model. The parameters are specified in the material file as given in Table 2.5.



**Model by Arora**

The model by Arora Arora _ is identified by ``mobility_formula`` = 2. It reads:

.. math::
   :label: dd_eq_dopdep2

      \mu = \mu_{min} +  \frac{\mu_d}{1+(N/N_0)^{A^*}} 
    
with

.. math::

    \mu_{min} = A_{min}(T/T_0)^{\alpha_m}, & \quad \mu_d = A_d(T/T_0)^{\alpha_d} \nonumber \\
    N_0 = A_N(T/T_0)^{\alpha_N}, & \quad A^* = A_a(T/T_0)^{\alpha_a} \nonumber 

The parameters are given in table at the end of the Chapter.



Field dependent mobility model
..............................


The field dependent mobility model describes the degradation of mobility at high driving
fields. It is identified by the identifier field_dependent. The electric field component
in direction of the current 
flow, the gradient of the electro-chemical potential or a so called "field parameter" can be
chosen as driving force:

  ``driving_force = efield | grad_fermi | field_parameter``

The default driving force is the gradient of the corresponding electro-chemical potential :math:`\nabla\phi` .
``field_parameter`` uses a field parameter given by :math:`\sqrt{E\cdot\nabla\phi}`  as driving force [Zakhleniuk]_ .

The model is based on the Caughey-Thomas model, refined by Canali [6]:

.. math::
   :label: dd_eq_fielddepmodel

    \mu = \frac{\mu_{lowfield}}{\left(1 + \left(\frac{\mu_{lowfield} |\mathbf{E}|}{v_{sat}}\right)^\beta \right)^{1/\beta}}
    
with

.. math::

    \beta = \beta_0(T/T_0)^b 

:math:`|E|` is the modulus of the driving field, multiplied by a damping factor :math:`n/(n+n_0)`, where :math:`n` is the electron or hole density, and :math:`n_0` is a parameterwith default :math:`10^{9}`.  :math:`\mu_{lowfield}` is the low-field mobility. For the latter
one can specify the model to be used using the parameter ``lowfield_model`` . As default
the doping dependent model is used.
There are two models for vsat, identified with ``Vsat_Formula = 1`` and 2. 

Formula 1 reads

.. math::
   :label: dd_eq_fielddepvel1

    v_{sat} = v_{sat,0} (T/T_0)^{-\gamma}
 

Formula 2 reads

.. math::
   :label: dd_eq_fielddepvel2

    v_{sat} = \max(A_{vsat} - B_{vsat} (T/T_0), v_{min})

The parameters for the field dependent mobility model are summarized in Table 2.7.


Field assisted mobility model
............................................


The field assisted mobility model describes the enhancement of the carrier mobility by an electric field in organic
semiconductors. It is identified by the identifier field_enhanced. 

The model is given by equation [devometterelareference]_:

.. math::
   :label: dd_eq_fieldassistedmodel

    \mu = \mu_0 e^{\sqrt{|E|/E_0}}

    
where :math:`|E|` is the modulus of the driving field, :math:`\mu_{0}` is the zero-field mobility and :math:`E_0` is a critical field strength.

The parameters for the field assisted mobility model are the following (summarized in Table :ref:`Field assisted mobility parameters<dd_field_assmob>`):

  ``mu0``
    The mobility at low electric field.

  ``E0``
    The critical electric field strength.


..  _dd_field_assmob :

.. math::
   :nowrap:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l|r|r|l}
    \hline
    \textbf{parameter name} & \multicolumn{2}{r|}{\textbf{default value}} & \textbf{units} \\
    \hline
    \hline
    \texttt{mu0} & \multicolumn{2}{r|}{0.0054} & cm$^2$/Vs \\
    \texttt{E0} & \multicolumn{2}{r|}{$3\cdot 10^5$} & V/cm \\
    \hline
    \end{tabular}
    \caption{Field assisted mobility parameters}
    \end{table}


Hopping mobility model
............................................

The hopping mobility model is obtained as a parametrical fit of numerical solutions of a system of master equations implementing the Miller-Abrahams hopping model between localized states with a gaussian energy distribution (see [Pasveer]_). It is identified by the identifier ``hopping_mobility``. The expression for mobility is:

.. math::

    \mu (T, \rho, F) = \mu(T,\rho) f(T,F)

.. math::

    \mu(T, \rho) = \mu_0 c_1 \exp \left[-c_2 \left(\frac{\sigma}{k_B T}\right)^2\right] \exp \left[\frac{1}{2 k_B T} \left(\frac{\sigma^2}{k_B T} - \sigma \right) \left(2 \rho a^3 \right)^\delta \right]

.. math::

    f (T,F) = \exp \left\{ 0.44 \left[\left( \frac{\sigma}{k_B T} \right)^{\frac{3}{2}} - 2.2 \right] \left[ \sqrt{1+0.8\left(\frac{eFa}{\sigma} \right)^2} -1 \right] \right\}

where :math:`\rho` is the carrier density (electrons or holes); :math:`F` is the electric field; :math:`T` the temperature :math:`a = N_0^{-\frac{1}{3}}` is the average distance between sites (:math:`N_0` is the maximum carrier density (:math:`cm^{-3}`) for the gaussian DOS); :math:`\sigma` is the variance of the gaussian DOS; :math:`\nu_0` is the attempt to jump frequency of the Miller-Abrahams model; :math:`c_1 = 1.8 \times 10^{-9}`; :math:`c_2 = 0.42` and

.. math::

   \delta = 2 \frac{\log\left(s^2 -s\right) - \log \left(\log 4\right)}{s^2} \ \ \ \ \mu_0 = \frac{a^2 \nu_0 e}{\sigma}

where :math:`s=\frac{\sigma}{k_B T}`. 

It is important to notice that this formula makes sense physically only if the condition :math:`\left( \frac{\sigma}{k_B T} \right)^{\frac{3}{2}} - 2.2 > 0` is fulfilled. In the TiberCAD input file the hopping mobility model has to be set providing values for :math:`\sigma`, :math:`N_0` and :math:`\nu_0` using the respective keywords ``sigma = [value]``, ``N0 =`` ... and ``nu0 =`` ... . This mobility model makes sense only if used together with the gaussian density of states, hence the values of ``sigma`` and ``N0`` have to be the same as the corresponding values set in the density of states model for the gaussian DOS.

As explained in [Santoni]_ two cut-off have been implemented. For :math:`F > \frac{2 \sigma}{ea}`

.. math::

   f(T,F) = f(T, \frac{2 \sigma}{ea})

For :math:`\rho > 0.1 N_0`

.. math::

   \mu (T, \rho) = \mu (T, 0.1 N_0)


Polarization models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For simulations involving materials with nonzero electric polarization (such as nitrides)
it is important to include the effect of polarization. This is done by specifying the models
for spontaneous (pyro-) and piezoelectric polarization using the keyword ``polarization``
with the types ``pyro`` and ``piezo`` ::

  polarization pyro {}
  polarization piezo {}

As for all models, if they do not have individual options, they can be specified together
by writing ``polarization (pyro, piezo) {}``.



Spontaneous polarization
............................................



The spontaneous polarization model (sometimes also called 'pyroelectric polarization') imposes a constant electric polarization P along
the symmetry-braking direction of the crystal. Crystals with wurtzite structure like Nitrides
have strong polarization fields along the c-direction. 
The value of the polarization usually is taken from the database, but it can be overridden from the input file by specifying
the option Pz, meaning the value of the spontaneous polarization along c-direction ([0001]). 
Alternatively, one can specify explicitly a polarization vector using the option
``P = (Px, Py, Pz)`` . This is useful to impose an arbitrary constant polarization field.



Piezopolarization
............................................



The piezoelectric polarization is strain induced and given by the linear relationship

.. math::
   :label: dd_eq_piezo

    P^{pz} = e_{ikl}\varepsilon_{kl}

where :math:`\varepsilon_{kl}` is the strain tensor. The piezoelectric coefficients :math:`e_{ikl}` are stored in the database.
The strain is obtained from the simulation specified in the ``Physics`` section, but it can
be overridden by providing a name for the strain simulation inside the polarization block
using the ``strain_simulation`` option. ::

  polarization (piezo) {strain_simulation = strain}



.. _DD_trapmodels:

Trap models
^^^^^^^^^^^^^^^^^^^^^^^^

Currently single level traps are implemented in TiberCAD. Traps can be normally neutral
or normally charged electron or hole traps, or a fixed charge. Common options for all
models are

  ``type`` 
         The type of traps. One of ``eNeutral``, ``hNeutral``,
         ``donor``, ``acceptor`` or ``fixed_charge``.
         (Only necessary if not provided as second keyword). 

  ``Nt`` 
         The trap density in |cm3| (or |cm2| for surface traps).

  ``Et`` 
         The trap level in eV with respect to the reference energy.

  ``reference`` 
         The reference energy. The default is ``m`` for midgap.
         Possible values are ``cb``, ``vb`` or ``m`` 

  ``recombination_center``
         Flag to switch on and off the recombination through a trap


For ``reference = m`` for example, the trap energy is given as :math:`E_{trap} = E_{midgap} + Et`. 
In the other cases it is :math:`E_{trap} = E_c - Et` or :math:`E_{trap} = E_v + Et`.
The following trap types are implemented:

  ``eNeutral``
      The trapped electron density is given by

      .. math::
       :label: dd_eq_eneutral

        n_t = \frac{N_t}{1 + \exp(\frac{E_{trap} - E_{F,n}}{k_BT})}


  ``hNeutral``
      The trapped hole density is given by

      .. math::
       :label: dd_eq_hneutral

        p_t = \frac{N_t}{1 + \exp(-\frac{E_{trap} - E_{F,p}}{k_BT})}


  ``donor``
      The density of ionized traps is given by

      .. math::
       :label: dd_eq_donor

       N^+_t = N_t - \frac{N_t}{1 + \exp(\frac{E_{trap} - E_{F,n}}{k_BT})}


  ``acceptor``
       The density of ionized traps is given by

       .. math::
        :label: dd_eq_acceptor
   
        N^-_t = N_t - \frac{N_t}{1 + \exp(-\frac{E_{trap} - E_{F,p}}{k_BT})}

If traps are specified, the total charge density in the Poisson equation is modified to
include the charged trap densities:

.. math::
   :label: dd_eq_totchargedensity

   \rho = e\left(p - n + N^+_D - N^-_A - \sum n_t + \sum p_t + \sum N^+_t - \sum N^-_t \right)

Additionally, each trap can induces a SRH recombination term of the form

.. math::
   :label: dd_eq_trapsrh

   R_t = N_t \frac{v_{th}^n\sigma^nv_{th}^p\sigma^p(np -n_i^2)}{v_{th}^n\sigma^n(n+n_1) + v_{th}^p\sigma^p(p+p1)}

where :math:`\sigma^{n,p}` are the capture cross sections, :math:`v_{th}^{n,p}` the thermal velocities and (for Boltzmann statistics)

.. math::
   :label: dd_eq_n1p1

   n_1 = n_{i,\mathrm{eff}}\exp(E_{trap}/k_BT),\quad p_1 = n_{i,\mathrm{eff}}\exp(-E_{trap}/k_BT)

The SRH (see :ref:`DD_recombinationmodels`) recombination model associated with a trap has to be enabled explicitly by using the option
``recombination_center = true`` in the trap definition.
If ``recombination_center = true`` is specified, trap-assisted tunneling (TAT)as described in the following section can be enabled by adding a subblock
``trap_assisted_tunneling`` as follows::

  trap_assisted_tunneling
  {
    tunneling_mass = 0.45
  }

   
This will activate the Hurkx TAT model with the specified tunneling mass, using the trap energy level given in the trap description.


Tunneling models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several simplified tunneling models are implemented in tiberCAD. They are derived from quantum mechanical
models under simplifying assumptions and cast into local recombination-generation models.




Trap-assisted tunneling
........................

Trap-assisted tunneling (TAT) in  both  forward and  reverse  bias is  taken  into  account with  a  recombination  model  proposed  by  Hurkx (see  [Hurkx]_). In this  model, a modified  expression  for  the  SRH recombination is  found,  where the carrier SRH lifetimes are modified by field-effect  functions :math:`\Gamma_n` and :math:`\Gamma_p` as :math:`\tau_{n,p} = {\tau_{n,p}}^0/(1+\Gamma_{n,p})`. The field-effect functions vanish  for weak  electric fields, yielding the conventional  SRH  formula.


The model for  the  Trap-assisted tunneling can  be enabled by defining a ``recombination`` submodel of  type ``srh``, adding the  keyword ``trap_assisted_tunneling`` as  follows::

  recombination srh 
  {
       Et = 0.56  
       reference = cb
       trap_assisted_tunneling = true
  }



The parameters  are  the  following (see also  section :ref:`DD_trapmodels`   )  : 

  ``Et``
      The trap level in eV with respect to the reference energy.

  ``reference`` 
         The reference energy. The default is ``m`` for midgap.
         Possible values are ``cb``, ``vb`` or ``m`` 

  ``trap_assisted_tunneling``
      if **true** tunneling  through the  defined trap is switched on 



(see [Hurkx]_).



Band-to-band tunneling
........................

For direct band-to-band tunneling, a model that is also due to  Hurkx is implemented. This model can be written as a local  generation-recombination process:


.. math::
   :label: Hurkx_model

   G^{b2b} = B E^\sigma exp(E_0/E) 


The model is activated by defining a ``recombination`` submodel of  type ``band2band`` as  follows::


  recombination band2band
  {
    B = 4e14  
    E0 = 1.9e7  
    sigma = 2.5 
  }


The parameters for the band to band  tunneling  model are summarized in Table :ref:`Hurkx band to band tunneling parameters<b2b_tunnel>`, together with  their default value:




..  _b2b_tunnel :

.. math::
   :nowrap:

    \begin{table}[!ht]
    \center
    \begin{tabular}{l|r|r|l}
    \hline
    \textbf{parameter name} & \multicolumn{2}{r|}{\textbf{default value}} & \textbf{units} \\
    \hline
    \hline
    \texttt{B} & \multicolumn{2}{r|}{$4\cdot 10^{14}$} & 1/cm$^3$s \\
    \texttt{E0} & \multicolumn{2}{r|}{$1.9\cdot 10^7$} & V/cm \\
    \texttt{sigma} & \multicolumn{2}{r|}{$2.5$} & 1 \\
    \hline
    \end{tabular}
    \caption{Hurkx Band to band tunneling parameters}
    \end{table}



(see [Hurkx]_).




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

    -\nabla[\mu_n n (\nabla\phi_n + P_n \nabla T) ]  & =  v_n(n - n_0) \\
    -\nabla[\mu_p p (\nabla\phi_p + P_p \nabla T) ]  & =  -v_p(p - p0)

The options ``zero_field`` , ``zero_grad_fermi_e`` and ``zero_grad_fermi_h`` can be used,
which when set to ``true`` will impose zero normal electric field and zero normal gradient
of the electron and hole electro-chemical potential, respectively. The latter are special
cases of surface recombination velocities ( :math:`v_{rec} = 0` ).

Contacts are defined by blocks with keyword Contact, for example::

  Contact anode 
  {
   type = ohmic
   [regions = (anode1, anode2)]
   voltage = $Vd
  }
    
An area factor can be specified for contacts using the keyword ``area_factor`` . The
contact current will be multiplied by this factor.

For interfaces and surfaces, the same syntax can be used (optionally one can  use the keywords ``Interface`` or ``Boundary``), however they do usually not need to be defined explicitly.

Ohmic contact
............................................



The ohmic contact (identifier ``ohmic``) has no further parameters.


Schottky contact
............................................



A Schottky contact (identifier ``schottky``) has the additional parameter ``barrier`` , which
signifies the energy difference between the semiconductor band edge and the fermi energy
in the metal. As default, the barrier is taken with respect to the conduction band. By
specifying ``band = v`` the barrier can be imposed with respect to the valence band (p-
type contact). Alternatively, the metal work function can be defined using the keyword
``work_function`` or the keyword ``metal_fermilevel``. The latter is just the work function with
inverted sign.

.. note:: The  value given in  ``work_function`` or ``metal_fermilevel`` has to be
          aligned with the band energies given in the material files,
          *not* with that resulting from simulation.

The ``fixed_barrier`` controls the behaviour of the barrier height for strained materials. 
If it is set to true, the barrier willbe independent of strain (default behaviour). 
If it is set to ``false`` , the given barrier is used as barrier for the unstrained case 
and will depend on strain during simulation. If the metal work function is specified, the barrier 
will be strain dependent as default. Thermionic emission is by default switched on, but can 
be disabled by specifying ``thermionic_emission = false`` .

.. warning::  If a Schottky contact is touching different materials, one should specify the work
              function instead of the barrier.

``barrier_lowering = true`` activates the image charge lowering of the Schottky barrier (see [Santoni]_);
if it is not specified, by default it is set to false. The Schottky barrier :math:`\phi_B` is lowered by a
quantity depending on the electric field :math:`F`:

.. math::
   \phi_B = \phi_{B0} - \sqrt{\frac{eF}{4 \pi \varepsilon}}

Internally TiberCAD sets conduction and valence band edges, then it fixes the barrier setting a value of the
fermi level inside the band gap. It is the important to notice that it is not possible to lower the barrier for
both electrons and holes, in fact lowering the barrier for electrons (``band = c``) means to shift the fermi level 
toward the conduction band, and accordingly the barrier for holes raises, and vice versa lowering the barrier for
holes (``band = v``) means to shift the fermi level toward the valence band, and accordingly the barrier for
electrons raises.

``scott_malliaras = true`` tells TiberCAD to use the [Scott]_ and Malliaras model of recombination velocity at the contact;
if it is not specified, by default it is set to false and the thermal velocity is used instead.

.. math::
   v_{rec} = v_{SM} \left(f\right) = \frac{v_{SM} \left(0\right)}{4} \left(\frac{1}{\psi^2 \left(f\right)} - f\right)

.. math::
   v_{SM} \left(0\right) = \frac{16 \pi \varepsilon \left(k_B T\right)^2 \mu}{e^3}

.. math::
   \psi \left(f\right) = f^{-1} + f^{-\frac{1}{2}} - f^{-1} \left(1 + 2 f^{\frac{1}{2}}\right)^{\frac{1}{2}}

.. math::
   f = \frac{eFr_C}{k_B T} \ \ \ \ r_C = \frac{e^2}{4 \pi \varepsilon k_B T}

The Scott and Malliaras model takes into account the effect of the electric field, and it is used in organic semiconductors
because it does not depend on the effective mass, a parameter not well defined in organic materials.


Interface/surface model
............................................




The free surface or interface model (identifier interface) can include surface charges
due to traps and surface recombination. Their definition can be found in section :ref:`DD_trapmodels`.

Each trap model will induce automatically a SRH recombination model as in the bulk
case.

:math:`Schr\ddot{o}dinger`/Poisson/Drift-Diffusion calculations
-----------------------------------------------------

tiberCAD is able to perform selfconsistent :math:`Schr\ddot{o}dinger`-Poisson or :math:`Schr\ddot{o}dinger`-Drift-Diffusion calculations.  For this purpose, a **density_of_states** model of type *quantum*
has to be specified in **Physics** for at least one of
the carriers,  and a *selfconsistent* simulation should be defined in the *Selfconsistent*
block (see :ref:`Selfcons`).  The following option, to be specified in the Physics section,  controls the behaviour of the ``selfconsistent`` simulation.

 ``use_density_predictor`` 
    When set to true, a predictor-corrector scheme will
    be adopted in the selfconsistent cycle. The Poisson/Drift-Diffusion solver does not
    just take the particle densities as given by the :math:`Schr\ddot{o}dinger` calculation, but it will
    assume a dependency of the density on the potentials of the form

    .. math::
     :nowrap:
     :label:

     \begin{equation}
     \rho(\varphi, \phi_n, \phi_p) = \frac{\rho_{\mathrm{quantum}}(\varphi^0, \phi_n^0, \phi_p^0)}{\rho_{\mathrm{classical}}(\varphi^0, \phi_n^0, \phi_p^0)}\rho_{\mathrm{classical}}(\varphi, \phi_n, \phi_p)
     \end{equation}

    where :math:`(\varphi^0, \phi_n^0, \phi_p^0)` are the potentials for which the quantum density was calculated.
    ``use_density_predictor = true`` is the preferred method for selfconsistent
    :math:`Schr\ddot{o}dinger`-Poisson/Drift-Diffusion calculations and is enabled by default.


For example, in **Module** *driftdiffusion* 
::

  Physics
  {    
   conduction_band 
    {
      density_of_states quantum
      { 
        # where to get the quantum density from
        quantum_simulation = quantum_el
        barrier_regions = buffer_quantum

	classical_DOS bulk_kp 
        {
          strain_simulation = strain
        }

      }
    }

  valence_band 
    {
      density_of_states quantum
      { 
        # where to get the quantum density from
        quantum_simulation = quantum_hl
        barrier_regions = buffer_quantum

	classical_DOS bulk_kp 
        {
          strain_simulation = strain
        }

      }
    }

In this example, we include both electrons and holes in the selfconsistent simulation, by defining a *quantum* model of DOS for both *conduction_band* and *valence_band*. Here ``quantum_el`` and  ``quantum_hl`` are  the  *efaschroedinger* simulations which  calculate the quantum densities for  electron and  holes, respectively.

Then, in **Module** *selfconsistent* ::

  Module selfconsistent
  {
   solve = (quantum_el, quantum_hl, driftdiffusion)
   # we do not use relaxation, but a predictor-corrector scheme
   #relaxation_factor = 0.5
   max_iterations = 10
   absolute_tolerance = 1e-3
   relative_tolerance = 1e-8
   monitor = true
   #xmonitor = true
  }



Again, ``quantum_el`` and  ``quantum_hl`` are  the  *efaschroedinger* simulations which  calculate the quantum densities for  electron and  holes. The  ``solve`` statement here specifies the order of execution in  the  self-consistent  cycle, which is  repeated until the  requested tolerance is  reached.




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



Example 1: pn diode
--------------

The following example shows a minimal Drift-Diffusion module definition for a pn junction.

----

::

  Module driftdiffusion
  {
    name = dd
    #regions = (pside, nside)
    plot = (Ec, Ev, eDensity, hDensity)

    Physics
    {
      recombination srh {}
        
      mobility doping_dependent {}
    }
        
    Contact anode 
    { 
      voltage = $Vd 
    }

    Contact cathode 
    { 
      voltage = 0 
    }
  }

----


.. _DD_Ex2:

Example 2: Mosfet
---------------

In this second example we show a 2D simulation of a silicon Mosfet device.
The  ``GMSH`` model (see :ref:`GMSH_Ex2`) is shown in Fig. :ref:`GMSH model of the Mosfet <fig_dd_mosfet>`.
 
.. _fig_dd_mosfet:

.. figure:: ../data/DDMosfetMesh.png
    :align: center
    :scale: 40%

    ``GMSH`` model of the Mosfet showing the mesh and the region labels

The model consists of a p-doped Si substrate (``substrate``), two highly n-doped access regions (``contact``),
a thin gate oxide (``oxide``) and source, gate, drain and back-side contacts.

We want to simulate a set of output characteristics and a transcharacteristic for this Mosfet, using two distinct input files: ``outputchar.tib`` and ``transchar.tib``.
In these two files we define only the ``Module sweep`` and ``Simulation`` blocks.
The device and model definitions are put into a third file ``mosfet.tib``, which is included in the other two files using the syntax

::

   @include mosfet.tib

The device definition found in ``mosfet.tib`` is shown in the following listing:

::

  Device mosfet
  {
    meshfile = mosfet.msh
  
    material = Si

    Region substrate 
    {
      Doping
      {
        density = 1e18
        type = acceptor
      }
    }
   
    Region contact
    {
      Doping
      {
        density = 5e19
        type = donor
      }
    }

    Region oxide
    {
      material = SiO2
    }
  }


The material is defined globally in the ``Device`` section.
For the oxide it has to be overridden in the correspondent ``Region`` block.

The followng shows the module definition for the Drift-Diffusion simulation:

::

  Module driftdiffusion
  { 
    #coupling = electrons

    plot = (Ec, Ev, eQFermi, eDensity, eCurrentDensity, eMobility,
            hQFermi, hDensity, hCurrentDensity, hMobility,
            NetRecombination, ElField, ElPotential, ContactCurrents)

    Solver
    {
      type = linesearch

      linear_solver
      {
        method = pconly
        preconditioner = lu
      }
    }

    Physics
    {
      
      recombination srh {}

      mobility
      {
        type = field_dependent
        low_field_model = doping_dependent
      }
    }

    Contact gate
    {
      type = schottky
      barrier_height = 3.0

      voltage = $Vg

      area_factor = 0.1
    }

    Contact source
    {
      voltage = 0.0
      area_factor = 0.1
    }
 
    Contact backcontact
    {
      voltage = 0.0
      area_factor = 0.1
    }

    Contact drain
    {
      voltage = $Vd
      area_factor = 0.1
    }
  }


The commented option ``coupling = electrons`` shows how a unipolar simulation can be set up.
This example however will be  simulated in bipolar mode.
The ``Solver`` defines options for the nonlinear solver (the shown options are the default ones).
In this case, a linesearch approach is used to refine the nonlinear Newton steps.
The linear solver for each Newton step (defined in the ``linear_solver`` block) uses a complete LU factorisation as preconditioner (``preconditioner = lu``).
In this case, the linear solver method can be chosen to be the application of the preconditioner only (``method = pconly``), instead of using an iterative
approach.

The ``Physics`` block contains the definition of a few physical models to be used.
For the particle density we use Fermi-Dirac statistics, which is the default.
We use Shockley-Read-Hall recombination since we solve for both electrons and holes.
For the mobility, we use a field-dependent model instead of the default constant mobility model.
The low-field mobility is chosen to be calculated from a doping dependent model (which is the default).

The contacts are defined in the ``Contact`` blocks.
For the gate we specify ``schottky`` as type (the default is ohmic contact), providing a suitable barrier height.
The ``area_factor = 0.1`` indicates that we assume a transistor with 1 mm gate width.

Next, we create a file ``transchar.tib`` containing the definitions for the simulation of the transcharacteristic, as given in the following listing::

  @include mosfet.tib

  Module sweep
  {
    name = sweep_drain
    solve = driftdiffusion

    variable = $Vd
    start = 0.0
    stop = 1.0 
    steps = 5 
  }

  Module sweep
  {
    name = sweep_gate
    solve = driftdiffusion

    variable = $Vg
    start = -0.5
    stop = 1.5
    steps = 100 

    max_step = 0.1
  
    plot_data = true
  }

  Simulation
  {
    solve = (sweep_drain, sweep_gate)
    resultpath = output_transchar 
    output_format = vtk
  }

On the first row we include the device definition from ``mosfet.tib``.
Then we define two sweeps ``sweep_drain`` and ``sweep_gate``.
The first one will ramp the drain to 1.0 V in 5 steps, without plotting the results.
The second one will then perform a sweep on the gate voltage from -0.5 V to 1.5 V,
plotting the results after each step and produciong the transfer characteristic.
The option ``max_step = 0.1`` limits the maximum voltage step to 0.1 V.
This is useful, since the solver has first to reach the initial gate voltage of -0.5 V, starting from 0 V.
Using this option it will do this in steps of 0.1 V.
In the ``Simulation`` block we simply have to specify the two sweeps in the correct order.

Running the simulation with ::

  tibercad transchar.tib
 
will produce in particular a set of files ``*.vtu`` for each 
step of the sweep ``sweep_gate``, and the file ``sweep_gate_driftdiffusion.dat`` containing the voltage-current
characteristics for each contact, shown in Fig. :ref:`fig_dd_mosfet_transchar`.

 
.. _fig_dd_mosfet_transchar:

.. figure:: ../data/DDMosfetTranschar.png
    :align: center
    :scale: 80%

    Mosfet transcharacteristic


For the simulation of the output characteristics we create a file ``outputchar.tib`` with the following content::

  @include mosfet.tib

  Module sweep
  {
    name = sweep_drain
    solve = driftdiffusion

    variable = $Vd
    start = 0.0
    stop = 2.0 
    steps = 40
  
    plot_data = true
  }

  Module sweep
  {
    name = sweep_gate
    solve = sweep_drain

    variable = $Vg
    start = 0.0
    stop = 1.5
    steps = 6
  }

  Simulation
  {
    solve = sweep_gate
    resultpath = output_outputchar  
  
    output_format = vtk
  }



As before, we include the device definition using the ``@include`` statement.
Then we define a sweep on the drain voltage with name ``sweep_drain`` and a second sweep ``sweep_gate``
to sweep the gate voltage.
In the latter, we specify ``sweep_drain`` in the ``solve`` option, creating thus a nested sweep.
For each gate voltage, a sweep over the drain voltage will be performed.

Running the simulation will produce a file for each couple of values (``$Vg``, ``$Vd``), and a file containing the output
characteristic for each value of ``$Vg``.
The resulting set of output characteristics is shown in Fig. :ref:`fig_dd_mosfet_outchar`.

 
.. _fig_dd_mosfet_outchar:

.. figure:: ../data/DDMosfetOutchar.png
    :align: center
    :scale: 80%

    Mosfet output characteristics




.. |cm2| replace::  cm\ :sup:`2`

.. |cm3| replace::  cm\ :sup:`3`


.. rubric::  Footnotes


