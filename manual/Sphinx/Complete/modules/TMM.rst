..   <marker>

.. _TMMTheory:

TMM
=================================================






Theory
------------

TMM modules is a Finite Element solver for analysis propagation of the light in the multi-layer
strucrue, and calculating generation rate, Intensity and Electric Field inside the layers. 
Combination of TMM and Drift and Difution modules Could be used in the Electro-Optical simulation of the Solar cells. 
It can also be used to calculate Reflection, Transmision and absorbtion of the multi-layer structures.

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

Incident wave
---------------
Incident wave boundry condition should be used to determin incoming wave point. It should be note that this point shouldn't be at the middle of the geometry. 

  		Contact point1
  		{
  			type = incident_wave
  		}


reflectivity
---------------
this parameter indicate reflectivity at the other side of Incoming wave. 
back_reflectivity = 1 is equal to total reflection at the end of the geometry.

		back_reflectivity = 0

incident angle
---------------
This parameter determin angle of the incoming wave. The value should be less than 90[deg].

		incident_angle = 0

down_lambda & up_lambda
---------------
These two parameter will detrmin range of wavelength which will be used in the calculations.
Usauly it should be from 300 to 1000(nm)

  wavelength_lower_lim = 300
  wavelength_steps = 10
  wavelength_uper_lim = 1000



Incoherency
---------------
for Coherence layers this parameter should be zero. by setting this parameter to one, corresponding layer will be considered
as an Incoherent layer. This feature is usfull for simulation of thick layers such as Glass.

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



 



