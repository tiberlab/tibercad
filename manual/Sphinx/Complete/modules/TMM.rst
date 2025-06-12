..   <marker>

.. _TMMTheory:

TMM
===






Theory
------

TMM modules is a Finite Element solver for analysis propagation of the light in the multi-layer
strucrue, and calculating generation rate, Intensity and Electric Field inside the layers. Ligth source can be either external source
or internal emission fo the dipole.
Combination of TMM and Drift and Difution modules Could be used in the Electro-Optical simulation of the Solar cells. 
It can also be used to calculate Reflection, Transmision and absorption of the multi-layer structures.

Solution/Plot variables
-----------------------

The solution variables available for plotting and for interaction with other modules are
given in :ref:`Plotting variables<TMM_solutions>`.





Solver section
--------------------

The ``Solver`` section of the ``TMM`` module refers to a linear solver.
See section :ref:`Linear_solver`  for details on linear solver options.


Physics section
--------------------


In the following we will describe all the physical models. As mentioned in the Introduction in section :ref:`InputFileGetting`, submodels can be restricted to a subset of simulation regions.




Incident wave (External source)
-------------------------------

Incident wave boundry condition(external source) should be used to determine incoming wave point. It should be noted that this point shouldn't be at the middle of the geometry. 

  		Contact point1
  		{
  			type = incident_wave
  		}


Dipole Coordinate (Internal source)
-----------------------------------

With this parameter one can specify coordinate of the dipoles, the value should mesh based (similer to geometry).
it is possiable to solve for both external source and

	dipole_coordinate = 200

Or


	dipole_coordinate = ( 100 200 300)



Dipole Modes (Internal source)
------------------------------

Dipoles polarization can be 'TE' or 'TM', one can solve for both polarizations by typing 'TEM'.
Dipole orientation can be 'V' or 'H', it would be possiable to solve for both of these two modes by specifying 'VH'.

examples of TM mode with Horizental dirrection:

  polarization = TM
  orientation = H

Solving for TE and TM mode with both horizental and vertical dirrection:

  polarization = TEM
  orientation = VH

Dipole Power (Used only for Internal source)
---------------------------------------------

This parmeter specify the power of the Dipole. in the case of multiple dipole, by specifying the scaler dipole power, power of all
dipoles would be identical, but it is possiable to defing non uniform power.

Uniform emission:
	  dipole_power = 0.2
	  dipole_coordinate = ( 100 200 300)

Non-Uniform emission:
	  dipole_power = (0.2 0.3 0.4)
	  dipole_coordinate = ( 100 200 300)

Dipole Radial Wave Number (Used only for Internal source)
----------------------------------------------------------

Emission of the dipole angle can be studied by varying radial wave number of the dipole emission.  And instead of giving the exact value
of the radial wave number, one could specify the ratio of the radial wave number('x' component) with respect to wave number in the dipole layer
('z' componenet). in this way ratio of '0' equal to normal incident and ratio of '1' means 90[deg] emission angle.

Solving for single angle:
	
	  dipole_steps = 0
	  dipole_ratio = 0 // solving only for 0[deg] emission

Solving for a range of angles:  
	  dipole_steps = 100
	  dipole_ratio = 1 // solving from 0[deg] to 90[deg] emission angle with 100 steps

reflectivity ( Used only for External source)
---------------------------------------------
this parameter indicate reflectivity at the other side of Incoming wave. 
back_reflectivity = 1 is equal to total reflection at the end of the geometry.

		back_reflectivity = 0

incident angle ( Used only for External source)
-----------------------------------------------
This parameter determin angle of the incoming external wave. The value should be less than 90[deg].

		incident_angle = 0

Wave Length ( Used for External source as well as Internal source)
------------------------------------------------------------------
There are two ways to specify emission angle of the external and internal source.
First by giving exact value of the Wave length in form of scaler or vector in "nm" unit:

	wavelengths = 450

Or
	
	wavelengths = (300 400 500)
	
Secondly it is possiable to define a range of wavelengths and steps of the wave lengths:

  wavelength_lower_lim = 300
  wavelength_steps = 1
  wavelength_uper_lim = 950
  



Incoherency ( Used only for External source)
---------------
for Coherence layers this parameter should be zero. by setting this parameter to one, corresponding layer will be considered
as an Incoherent layer. This feature is usfull for simulation of thick layers such as Glass. Thisfeatre is only applied for external source
calculations.

  		Region Glass 
  		{
    			Incoherency= 1
    			material = Glass
  		}







Example 1
--------------
In the following example, we model a silver back contact and calculate optical parameters.
The model consits of 3 regions, one silver layer between two air layer and all layers are coherent.

	Device
	{
  		meshfile = bulk.msh 
  		mesh_units = 1e-9
  		dimension = 1
  		Region AIR1 
  		{
    			Incoherency = 0
    			material = Air
  		}
  		Region Ag
  		{
    			Incoherency= 0
    			material = Ag
  		}
  		Region AIR2 
  		{
    			Incoherency = 0
   			material = Air
  		}
	}

In the modules section, we will define which Tmm Outputs we are intested in for instance, Generation Rate  or Reflectio or... .

	Module tmm
	{ 
 		name = tmm_study 
  		plot = (GenerationRate,Intensity,ElectricField, Reflection,Transmission,Absorbtion)
	}

Then we have define wavelength spectrum for simulation by given value to two parameters of "up_lambda" and "down_lambda"(unit is nm).
make sure than "down_lambda" should has a lower value.

  		wavelength_lower_lim = 300
  		wavelength_steps = 10
  		wavelength_uper_lim = 1000

We should also define incomin wave point( in this example due to symetry it dosnt make any diffenet)

  		Contact point1
  		{
  			type = incident_wave
  		}

At the end we will define Normal Incident and zero back reflection.
  
		incident_angle = 0
		back_reflectivity = 0











Example 2
------------
In this example transmition and reflection of a thcik layer of Glass has been simulated using TMM module. The Glass is simulated by adding a Inoherence model for the layer. Incogerence parameter of the Glass is set to be '1' while other layer remain '0'.

	Device
	{
  		meshfile = bulk.msh 
  		mesh_units = 1e-9
  		dimension = 1
  		Region AIR1 
  		{
    			Incoherency = 0
    			material = Air
  		}
  		Region Glass 
  		{
    			Incoherency= 1
    			material = Glass
  		}
  		Region AIR2 
  		{
    			Incoherency = 0
   			material = Air
  		}
	}
in the modules section, wavelength spectrum has been set to be 300~950nm, incident angle and reflectivity are set to be '0'.
Point1 is the postion of the incoming wave.

	Module tmm
	{ 
 		name = tmm_study 
 		plot = (GenerationRate,Intensity,ElectricField, Reflection,Transmission,Absorbtion)
  		back_reflectivity = 0
  		incident_angle = 0
  		wavelength_uper_lim = 950
  		wavelength_lower_lim = 300
  		wavelength_steps = 10
 		Physics 
 		{

		}
 		Contact point1
  		{
 			type = incident_wave
 		}
	}











Example 3
------------
In this example a complete structure of a Perovskite solar cell is simulated using TMM module. 
Oerder of the layers are as follow:
Air / Glass / ITO / PTAA / PSK / SO2 / Ag / Air

The Device is defined as follow:
	Device
	{
  		meshfile = bulk.msh 
  		mesh_units = 1e-9
  		dimension = 1

  		Region AIR1 
  		{
  			Incoherency = 0
  			material = Air
  		}

  		Region Glass 
  		{
  			Incoherency = 0
  			material = Glass
  		}

  		Region ITO 
  		{
  			Incoherency= 0
  			material = ITO
  		}

  		Region PTAA 
  		{
  			Incoherency = 0
  			material = PTAA
  		}

  		Region PSK 
  		{
  			Incoherency = 0
  			material = PSK
  		}

  		Region SnO2 
  		{
  			Incoherency = 0
  			material = SnO2
  		}

  		Region Ag
  		{
  			Incoherency = 0
  			material = Ag
  		}

  		Region AIR2
  		{
  			Incoherency = 0
  			material = Air
  		}

	}

TMM Module configuration is:
	Module tmm
	{ 
 		name = tmm_study 
 		plot = (GenerationRate,Intensity,ElectricField, Reflection,Transmission,Absorbtion)
  		back_reflectivity = 0
  		incident_angle = 0
  		wavelength_uper_lim = 950
  		wavelength_lower_lim = 300
  		wavelength_steps = 10
 		Physics 
 		{

		}
 		Contact point1
  		{
 			type = incident_wave
 		}
	}



 Example 4
------------
This example solve for TE mode with H orientation for 0[deg] up to 90[deg] emission angle with dipole power of 0.2. Dipole located at coordinate 500.

Module tmm
{ 
  
  name = tmm_study 
  plot = (Internal_Source_ElectricField,Internal_Poynting,Internal_Power,Internal_Absorption)
  wavelengths = 300
  dipole_steps = 200
  dipole_ratio = 1
  polarization = 1
  orientation = 0
  dipole_power = 0.2
  dipole_coordinate = 500

  Physics 
  {

  }

}


