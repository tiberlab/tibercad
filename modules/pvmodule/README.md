# pvmodule

Tibercad module *pvmodule* simulates IV characteristics of a solar module based on a 2D description of its layout. It uses a discretization into a lumped element circuit, including the sheet restsitivities of top and bottom contact layers. The elementary cell can be given as equivalent circuit, or as numerical data. The latter is transformed into a piecewise linear source.

The obtained discrete circuit is written as a Spice netlist, and solved using ngspice.

The module includes a framework for defining degradation models that can use local quantities form other modules, especially *wateringress*.

**NOTE:** module under development.
