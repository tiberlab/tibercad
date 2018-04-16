
..  _installGetting:

############
Introduction
############


Installation instructions
===================================

In the following, ``installpath`` denotes the directory where tiberCAD  gets installed.
A recent version  of ``GMSH`` (http://www.geuz.org/gmsh) will be installed together with
tiberCAD.

Prerequisites
-------------------

Get the installer package for your OS/architecture from http://www.tiberlab.com or
by contacting support@tiberlab.com. Table :ref:`Installer packages<input_installerpacks>`  lists the packages available for download.
To run tiberCAD you will also need a license file that you will have to copy into the
installation directory of tiberCAD.

.. In the Windows version, some graphical features such as graphical convergence 
.. monitors are only available if an X Window server is installed and running.

Windows installation procedure
--------------------------------------

To install tiberCAD in Windows, please run the setup program. 
``tibercad-3.0.0-i686-w64-mingw32_setup.exe`` .

During the installation you can choose the installation directory. 
After finishing installation, copy your license file ``tibercad.lic`` into the ``license`` subdirectory of the tiberCAD
installation directory (``installpath/license``), without changing the filename.


.. _input_installerpacks :

..  math::
    :nowrap:

    \begin{table}[!h]
    \center
    \begin{tabular}{l|c}
    \multicolumn{2}{c}{\textbf{Installer}} \\
    \hline
    \textbf{installer package name} & \textbf{Target architecture} \\
    \hline
    \hline
     &    \\
    \texttt{tiberCAD-3.0.0-i686-w64-mingw32\_setup.exe} & Windows 32-bit     \\
    \texttt{tiberCAD-3.0.0-x86\_64-linux\_installer.sh} & Linux 64-bit \\
    \hline
    \end{tabular}
    \caption*{Installer Packages}
    \end{table}



..    \texttt{tiberCAD-3.0.0-i686-linux\_installer.sh} & Linux 32-bit \\


Linux installation procedure
------------------------------------

To install tiberCAD under Linux, download and run the self-extracting installer
``tibercad-3.x-x86_64-linux_installer.sh``

.. ``tibercad-3.x-ARCH_installer.sh``, where ARCH corresponds to your hardware architecture,
and follow the installation instructions.

After installation, copy your license file ``tibercad.lic`` into the ``license`` subdirectory
of the tiberCAD installation path (``installpath/license``) without changing the filename. 
You can also provide the license file during installation.

tiberCAD is launched by means of a shell script, which is
installed together with the tiberCAD executable. It takes care of setting all necessary 
environment variables. 

If you have to call the executable directly, you have to
set the environment variable ``TIBERCADROOT`` to the tiberCAD installation directory (``installpath``).

Quick start guide
----------------------

In the ``examples`` subdirectory you can find several examples ready to run. More examples may be
available on http://www.tiberlab.com or http://www.tibercad.org.

tiberCAD has the following command line options:

  -v     print the version number and exit
  -b     run in batch mode, without asking for user input.

The ``-b`` option is useful only for the windows version, which attends a keystroke by the user before
exiting.

Windows
^^^^^^^^^^^^^^^^^^^^^

Open Windows Explorer and go to the tiberCAD installation directory. If you have
write permission in the installation directory, you can browse to an example directory
and start the simulation by double clicking the input file, e.g. ``bulk.tib``.

.. in :ref:`tut0step1`  

If not, copy the whole directory to a location in your personal area and run the examples
from there.

If you cannot run tiberCAD by double clicking a tiberCAD input file (the files with extension ``.tib``), then the input
files may not be correctly associated with the tiberCAD executable. In this case,
try to establish the association by right-clicking the input file, choosing 

  ``open with... >> Choose Program... >> Browse...``,
browsing to the tiberCAD installation directory and choosing the tiberCAD executable, ``tibercad.exe``.

Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After a correct installation of tiberCAD you should be able to run tiberCAD from the
command line using the command ``tibercad``. If not, you probably have to add the ``bin``
subdirectory of the tiberCAD installation path to your ``PATH`` environment variable or
start the tiberCAD executable using the absolute path (``installpath/bin/tibercad``).

Copy the directory of the example you want to run to your home directory
or any place you have write permissions for. Change to the example directory and
run tiberCAD by e.g. 

  ``$ tibercad bulk.tib``

.. (assuming :ref:`tut0step1` )

 


Bug reports / Feedback
-------------------------------------

Please send bug reports, feedback or suggestions to support@tibercad.org. When
submitting bug reports, please always include the full version number of tiberCAD you
are running. The full version number appears in the first line of output when running
the program:

::

  $ tibercad -v
  tiberCAD version 3.0.0 (x86_64-linux)





.. _GMSHTutorialGetting:


GMSH   Quick Tutorial
=================================================

To  use  tiberCAD,   as a first step you  need  to model the device and  generate  a suitable mesh  grid. This  can  be  done  by  using DEVISE module of ISE-TCAD 9.5 software package or GMSH program.


In the following we will see how to write a basic GMSH script in 1 and 2D; for any details please refer to GMSH manual  (http://geuz.org/gmsh/).


.. _GMSH_Ex1:

GMSH Example 1D 
-------------------------------------------


Step 1: Modeling the device
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Here we  will  refer  to  the  Example 0 *Bulk Silicon  in  1D* in the  Example directory. 
See :ref:`Input_Ex1` for the  description of the corresponding input file.

In a GMSH script, several variables can be defined and given a value in this way::


  L = 1;
  d = 0.01;


these are valid GMSH variables: ``L`` is just the length of the Si sample; ``d`` is the value of a *characteristic mesh length* (see below).



**GMSH modeling strategy: bottom-up  design**

In   gmsh  the idea  is  to   design  the  model  with a "bottom-up" strategy.
So,   first, points are  defined, then lines  connecting  points,  surface  connecting  lines,  and  so  on. No  superimposing  objects are  allowed.
This  means  that, once  defined  your  points,  you  may  connect them  with  lines but
different lines  must  not  have  parts  in  common (just  points);
the  same  works for  surfaces:  they  may  have only  lines  in  common,   but  no  intersections between  surfaces  are allowed.

.. warning::
             If  a  geometrical  model  with  not  null  intersections  between entities (points, lines,  surfaces, volumes) is  created, unpredictable results may  occur (gmsh  crashes  during  meshing, a  mesh is  created  which is  not  valid, etc.). 





**Definition of the geometrical entities Points**

::

  Point(1) = {0, 0, 0, d};
  Point(2) = {L, 0, 0, d};


In the definition of a geometrical point, the  first three expressions inside the braces on the right hand side give the three X, Y and Z coordinates of the point; the last expression ``d`` sets the *characteristic mesh length* at that point, that is the *size* of a mesh element, 
defined as the length of the segment for a line mesh element, the radius of the circumscribed circle for a triangle mesh element and the radius of the circumscribed sphere for a tetrahedron mesh element.


Thus, the smaller is the value of ``d``, the greater is the mesh density close to that point. The size of the mesh elements will then be computed in GMSH by linearly interpolating these characteristic lengths in the whole mesh.

.. warning::
             In a 1D simulation it is  assumed that the geometrical model is  restricted to the ``x`` axis. 
             Any other geometrical orientation  could  give unpredictable results.



**Definition of a geometrical entity Line**


::  

  Line(1) = {1,2};

The two expressions inside the braces on the right hand side  give the identification numbers of the start and end points of the line.



**Definition of the physical entity Physical Line ``bulk`` and of two physical entities Physical Point** 





Convenient *Physical Names*  are  to  be  assigned to  the  Physical entities. *Physical Names* consist of  strings  enclosed between  quotation marks.
The  syntax is the  following:

::
  
  Physical Line("bulk") = {1}
  ............
  Physical Point("Anode") = {1};
  Physical Point("Cathode") = {2};





The expression(s) inside the braces on the right hand side  give the identification numbers of all the geometrical lines that need to be grouped inside the *Physical Line*  or *Physical Point* .
In this way, in general, *physical regions* are created which associate together geometrical regions, and then the related mesh elements, which share some common physical properties. It's only these physical regions which can be referred to outside GMSH. In tiberCAD, this is done by associating one or more physical regions to a tiberCAD region through the keywords *Region* and  *mesh_regions* (see :ref:`Input File<InputFileGetting>`).




.. warning::
             In general, in a n-Dimension (``nD``) simulation, ``(n-1)D`` physical regions
             (points in 1D, lines in 2D, surfaces in 3D) are used by tiberCAD to impose the required 
             boundary conditions. Each  ``(n-1)D``  physical region defined in this way in GMSH will be associated                          in tiberCAD to a boundary condition (Contact) region. 
             Thus, in this case, Physical points *Anode* and *Cathode* will be associated respectively to two *Contacts* 
             (see :ref:`Input_Ex1`).

            


 


Step 2: Meshing the device
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



The ``.geo`` script file with the geometrical description can be run in GMSH, to display the modelled device and to mesh it through the GMSH graphical interface.
To  generate  the  mesh,  select  ``Mesh`` in  the main menu  of  GMSH  and  click  on  ``1D``, ``2D`` or ``3D`` depending  on  the  dimension of  your  simulation.
This  will  create  a  file  .msh  in  your  working  directory.


Alternatively, a ``non-interactive`` mode is also available in GMSH, without graphical user interface. For example, to mesh this 1D tutorial in non-interactive mode, just type in the command line ::

  gmsh bulk.geo  -1 -o bulk.msh 

where ``bulk.geo``  is the geometrical description of the device with GMSH syntax;
``-1`` means 1D mesh generation;

some command line options are::

  -1, -2, -3 

to perform 1D, 2D or 3D mesh generation, respectively.

::

  -o  mesh_file.msh 

to specify the name of the mesh file to be generated

In this way, a ``.msh`` has been generated and is ready to be read in tiberCAD. 




.. _GMSH_Ex2:

GMSH Example 2D
-------------------------------------------


In this  second  example  we  will  refer  to  the  Example 4 that you can 
find in the Example directory  (*mosfet.geo*). See :ref:`Input_Ex2` for  a  description of  the  Input file.



Modeling the device
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Again, as a first step, we have to model the device. 

Geometrical *Points* and *Lines* are  defined to design the device  structure; the  fourth parameter in *Point* assignement is  the   **characteristic length** associated to that point: this  is  an  essential feature to control the  mesh density and refine it where  necessary (usually in the channel region).   

.. warning::
             In a 2D simulation it is  assumed that the geometrical model is  restricted to  the  
             ``xy-plane (z = 0)``. Any  other geometrical orientation  could  give impredictable results


 
::
  
  Point(1) = {0, -h, 0, lsub};
  Point(2) = {0, 0, 0, lc};
  Point(3) = {xmax,-h,0.0,lsub};
  Point(4) = {-xmax,-h,0.0,lsub};
  Point(5) = {xmax,0,0.0,lh};
  Point(6) = {-xmax,0,0.0,lh};
  ..........................
  Line(1) = {4,1};
  Line(2) = {3,13};
  Line(6) = {4,14};
  Line(7) = {10,9};
  Line(8) = {12,2};
  Line(9) = {8,7};
  Line(10) = {11,8};
  Line(11) = {9,12};
  Line(13) = {7,6};
  ..........................


**Definition of a surface**

First a *line loop* is composed, listing all the  lines constituting the  boundary of the surface; then this  line  loop is  assigned to a  *Plane Surface* object (this  procedure can be alternatively performed through the  graphical interface).                          

::


  Line Loop(40) = {28,2,-34,33,8,29,-31,-30,-6,1};
  Plane Surface(41) = {40};
  ..........................


The obtained geometrical  surface is shown n Fig. :ref:`Surface<geomsurf>`


..  _geomsurf :

..  figure:: ../data/geo_surf.png
    :align: center
    :scale: 140%

    Surface




**Definition of  the Physical Surfaces**

Each of the  *Physical Surfaces* is composed by one or more geometrical *Plane Surface*. For example, *Physical surface* **contact** comprises in one   single physical region the two separated contact geometrical regions, while *Physical surface* **oxide** corresponds to the  oxide  region. 
The  *Physical surfaces* are the 2D Physical regions of  the  mesh and will  be  assigned to the related tiberCAD regions through the keyword *Region* and *mesh_regions*. (See :ref:`Input_Ex2`) ::

  Physical Surface("substrate") = {41}; // n-Si
  Physical Surface("contact") = {44,47}; // n+-Si
  Physical Surface("oxide") = {46}; // SiO2







                                               
**Definition of the Phisical Lines**

In this 2D simulation, 1D physical regions are used to carry information about boundary condition regions. In  other words, each *Phisical Line* corresponds to a boundary condition (a contact in the case of a driftdiffusion calculation). Thus *Physical Line* **source** refers to the source contact, *Physical Line* **gate**  to the gate contact, *Physical Line*  **drain**  to the drain contact.
The names of these *Phisical Lines*  will be  assigned to tiberCAD *Contacts*.

::

  Physical Line("source") = {13}; // source
  Physical Line("gate") = {39,38}; // gate
  Physical Line("drain") = {19}; // drain


The final geometrical model  is shown in Fig. :ref:`Geometrical model<geomodel>`


..  _geomodel :

..  figure:: ../data/geomosfet.png
    :align: center
    :scale: 140%

    Geometrical model 




Meshing the device
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``.geo`` script file with the geometrical description can be run in GMSH, to display the modelled device and to mesh it through the GMSH graphical interface.
Alternatively, a textual mode is also available in GMSH, without graphical user interface. For example, to mesh this 2D tutorial in non-interactive mode, just type:

::

  gmsh mosfet.geo  -2 -o mosfet.msh 


The final meshed model  is shown in Fig. :ref:`2D meshing<mesh>`


..  _mesh :

..  figure:: ../data/meshmosfet.png
    :align: center
    :scale: 140%

    Meshed model 


