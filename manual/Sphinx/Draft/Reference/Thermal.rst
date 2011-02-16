.. _Thermal:

Thermal
=================================================

Heat Source/constant
Kind: bulk 
Name: HeatSource
Type: constant
database_section: none
Options: 
H (double,0.0)

Heat Source/joule
Kind: bulk 
Name: HeatSource
Type: joule
database_section: none
Options: 
dd_simulation (string,”driftdiffusion”)

Boundary/Heat reservoir
Kind: interface
Name: Contact
Type: heat_reservoir
database_section: none
Options: 
temperature (positive double,300)

Boundary/Surface Thermal Resistance
Kind: interface
Name: Contact
Type: surface_resistance
database_section: none
Options: 
r_surf (positive double 0.0)
temperature (positive double (K),300)

Boundary/Flux
Kind: interface
Name: Contact
Type: flux
database_section: none
Options: 
   Flux  (double vector (W/cm2, (0.0,0.0,0.0))
