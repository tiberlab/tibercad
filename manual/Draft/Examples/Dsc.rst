.. _Dsc:


Simulation Dye Solar Cells
==============================

Model DSC
^^^^^^^^^

The DSC model is tagged as *dssc* . In **options** subsection we indicate the simulation name::

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
    

Generation model
^^^^^^^^^^^^^^^^

::

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
    

the coming of light and direction of rays.

Device
------

::

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
    

Physics
-----------

::

  dssc
    {
     generation = dssc_generation
    }
    

Solver
----------

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
    