tibercad
========

tiberCAD is a multiscale/multiphysics TCAD software, developed at the Electronic Engineering
Department of Tor Vergata University of Rome since 2005. Its main distinctive feature is the
combination of continuum models with atomistic ones in a single simulation environment.
The software is divided into a core and modules. The modules implement the actual physical models.
The code is heavily based on the [libMesh FEM library](https://libmesh.github.io), and uses [PETSc](https://www.mcs.anl.gov/petsc) and [SLEPc](http://slepc.upv.es/) numerical libraries.

See CONTRIBUTORS for a list of peope that contributed to tiberCAD.

License
-------

The tiberCAD core is distributed as is, under the conditions of the Lesser GNU General Public
License (LGPL) version 3.0.

Modules are released under the GNU General Public License (GPL) version 3.0, unless otherwise noted


Compilation
-----------

tiberCAD can be compiled by defining a site configuration in compilation/site_config.sh and running build_all. This will download necessary libraries and configure and compile all code. Example site configurations are in the site_configuration subdirectory.
Note that currently standard installations of PETSc and SLEPc, and therefore libmesh, cannot yet be used. In future, the compilation scripts might be changed to make this possible.

Tagging and releases
--------------------

The script tools/release.sh can be used to create release versions consistent with tags. Release tags are of the form v3.5.1, the script takes the version numer without the leading "v" as input.
