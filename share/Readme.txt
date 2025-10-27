
tiberCAD Multiscale Device Simulator
------------------------------------

Last updated: $Date$
Current Version: <TIBERVERSION>

Web   : http://www.tibercad.org
E-Mail: support@tibercad.org




Contents:

 1 Installation instructions
   1.1 Prerequisites
   1.2 Windows
   1.3 Linux
    
 2 Quick start guide
   2.1 Windows
   1.2 Linux

 3 Bug reports / Feedback



1 Installation instructions

In the following, VERSION denotes the version number of the TiberCAD release
you downloaded and INSTALLPATH denotes the directory where TiberCAD gets installed.
Version 2.5.0 of GMSH (http://www.geuz.org/gmsh) will be installed together with
TiberCAD. 


1.1 Prerequisites

Get the installer package for your OS/architecture from http://www.tiberlab.com
or by contacting support@tibercad.org. To run TiberCAD you will also need a license
file that you will have to copy into the installation directory of TiberCAD.


1.2 Windows

To install tiberCAD in Windows, run the setup program `tibercad-2.0.0_setup.exe'.
During the installation you can choose the installation directory. After finishing
installation, copy your license file (tibercad.lic) into the `license' subdirectory
in the TiberCAD installation directory (INSTALLPATH/license), without changing
its filename.


1.3 Linux

To install tiberCAD under Linux, download and run the self-extracting installer
`tibercad-2.0.0-ARCH_installer.bin', where ARCH corresponds to your machine architecture, 
and follow the installation instructions.

After installation, copy your license file `tibercad.lic' into the `license' subdirectory
of the tiberCAD installation path (`installpath/license') without changing the filename. 
You can also provide the license file during installation.

tiberCAD is launched by means of a shell script, which is
installed together with the tiberCAD executable. It takes care of setting all necessary 
environment variables. 

If you have to call the executable directly, you have to
set the environment variable `TIBERCADROOT' to the tiberCAD installation directory (`installpath').
 



2 Quick start guide

In the `examples' subdirectory you can find several examples ready to run. More examples may be
available on http://www.tiberlab.com or http://www.tibercad.org.

tiberCAD has the following command line options:

  -v     print the version number and exit
  -b     run in batch mode, without asking for user input.

The ``-b`` option is useful only for the windows version, which attends a keystroke by the user before
exiting.


2.1 Windows

Open Windows Explorer and go to the tiberCAD installation directory. If you have
write permission in the installation directory, you can browse to an example directory
and start the simulation by double clicking the input file, e.g. bulk.tib in Tutorial 1. If
not, copy the whole directory to a location in your personal area and run the examples
from there.

If you cannot run tiberCAD by double clicking a tiberCAD input file (the files with extension .tib), then the input
files may not be correctly associated with the tiberCAD executable. In this case,
try to establish the association by right-clicking the input file, choosing 

  open with... >> Choose Program... >> Browse...,

browsing to the tiberCAD installation directory and choosing the tiberCAD executable tibercad.exe.


2.2 Linux

After a correct installation of tiberCAD you should be able to run tiberCAD from the
command line using the command `tibercad'. If not, you probably have to add the `bin'
subdirectory of the tiberCAD installation path to your PATH environment variable or
start the tiberCAD executable using the absolute path (`installpath/bin/tibercad').

Copy the directory of the example you want to run to your home directory
or any place you have write permissions for. Change to the example directory and
run tiberCAD by (assuming Tutorial 1)

  $ tibercad bulk.tib




3 Bug reports / Feedback

Please send bug reports, feedback or suggestions to support@tiberlab.com. When
submitting bug reports, please always include the full version number of tiberCAD you
are running. The full version number appears in the first line of output when running
the program:

  $ tibercad -v
  TiberCAD version 2.0.0 (x86_64-linux)




