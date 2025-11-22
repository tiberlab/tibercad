tibercad
========

tibercad is a multiscale/multiphysics TCAD software, developed at the Electronic Engineering
Department of Tor Vergata University of Rome since 2005. Its main distinctive feature is the
combination of continuum models with atomistic ones in a single simulation environment.
The software is divided into a core and modules. The modules implement the actual physical models.
The code is heavily based on the [libMesh FEM library](https://libmesh.github.io), and uses [PETSc](https://www.mcs.anl.gov/petsc) and [SLEPc](http://slepc.upv.es/) numerical libraries.


License
-------

The tibercad core is distributed as is, under the conditions of the Lesser GNU General Public
License (LGPL) version 3.0.

Modules are released under the GNU General Public License (GPL) version 3.0, unless otherwise noted


Compilation
-----------

tibercad can be compiled by defining a site configuration in compilation/ and running build_all. This will download necessary libraries and configure and compile all code.
