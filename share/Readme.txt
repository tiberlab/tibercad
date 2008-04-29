
TiberCAD Multiscale Device Simulator
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
     1.3.1 Debian
     1.3.2 Other
    
 2 Quick start guide
   2.1 Windows
   1.2 Linux

 3 Bug reports / Feedback



1 Installation instructions

In the following, VERSION denotes the version number of the TiberCAD release
you downloaded and INSTALLPATH denotes the directory where TiberCAD gets installed.


1.1 Prerequisites

Get the installer package for your OS/architecture from http://www.tibercad.org
or by contacting support@tibercad.org. To run TiberCAD you will also need a license
file that you will have to copy into the installation directory of TiberCAD.


1.2 Windows

To install TiberCAD in Windows, run the setup program `tibercad_setup.exe'.
During the installation you can choose the installation directory. After finishing
installation, copy your license file (tibercad.lic) into the `license' subdirectory
in the TiberCAD installation directory (INSTALLPATH/license), without changing
its filename.

  You should now be ready to run your first TiberCAD example.


1.3 Linux

The installation procedure for the Linux version of TiberCAD depends on your
distribution. Download the installer package that best fits your setup.

1.3.1 Debian

A package for Debian 4.0 `etch' for the i386 architecture is provided (.deb).
Install it as root using dpkg or similar:

  # dpkg --install tibercad_VERSION_ARCH.deb

  The package will be installed in /usr/share/tibercad-VERSION. Next, copy your
license file (tibercad.lic) into /usr/share/tibercad-VERSION/license/ without
changing the filename.

  You should now be ready to run your first TiberCAD example.

NOTE: The debian version of TiberCAD depends on the following Debian packages:
   - libboost-regex1.33.1
   - libboost-filesystem1.33.1
   - libblas.so.3 (provided e.g. by atlas3-base)
   - liblapack.so.3 (provided e.g. by atlas3-base)

1.3.2 Other

If you have a distribution other than Debian 4.0 or you want to install TiberCAD
into a different directory, then use the .tgz or .tbz installation packages.
Unpack the archive, cd to the unpacked directory tibercad-VERSION and run the
install script.

  After installation, copy your license file (tibercad.lic) into the `license'
subdirectory of the TiberCAD installation directory (INSTALLPATH/license) without
changing the filename.

  You should now be ready to run your first TiberCAD example.

  

2 Quick start guide

2.1 Windows

Open Windows Explorer and go to the TiberCAD installation. In the `examples'
subdirectory you can find several examples ready to run. They are the same as
the tutorials on http://www.tibercad.org/documentation/tutorial/list.

  If you have write permission in the installation directory, you can go to
e.g. bulk_Si and start the simulation by double clicking the bulk.tib input file.
If not, copy the whole directory to a location in your personal area and run the
examples from there.

  If you cannot run TiberCAD by double clicking an input file (*.tib), then the
input files are probably not correctly associated with the TiberCAD executable.
In this case try to establish the association by right-clicking the input file,
choosing `open with... -> Choose Program... -> Browse...' and choosing the
TiberCAD executable.


1.2 Linux

In the INSTALLDIR/examples directory of the TiberCAD installation you can find
several examples ready to run. They are the same as the tutorials on
http://www.tibercad.org/documentation/tutorial/list.

  Copy the directory of the example you want to run, e.g. bulk_Si, to your home
  directory or any place you have write permissions for. Change to the newly
  created directory and run TiberCAD by (assuming the bulk_Si example)

  $ tibercad bulk.tib




3 Bug reports / Feedback

Please send bug reports, feedback or suggestions to support@tibercad.org.
