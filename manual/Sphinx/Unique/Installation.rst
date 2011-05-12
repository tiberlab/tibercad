
..  _installGetting:


Installation instructions
===================================

In the following, ``installpath`` denotes the directory where tiberCAD 2.0.0 gets installed.
Version 2.5.0 of ``GMSH`` (http://www.geuz.org/gmsh) will be installed together with
tiberCAD.

Prerequisites
-------------------

Get the installer package for your OS/architecture from http://www.tiberlab.com or
by contacting support@tibercad.org. Table 1 lists the packages available for download.
To run tiberCAD you will also need a license file that you will have to copy into the
installation directory of tiberCAD.

.. In the Windows version, some graphical features such as graphical convergence 
.. monitors are only available if an X Window server is installed and running.

Windows installation procedure
--------------------------------------

To install tiberCAD in Windows, run the setup program ``tibercad-2.0.0_setup.exe`` .

During the installation you can choose the installation directory. 
After finishing installation, copy your license file ``tibercad.lic`` into the ``license`` subdirectory of the tiberCAD
installation directory (``installpath/license``), without changing the filename.



..  math::

    \begin{tabular}{l|c}
    \multicolumn{2}{c}{\textbf{Installer}} \\
    \hline
    \textbf{installer package name} & \textbf{Target architecture} \\
    \hline
    \hline
     &    \\
    \texttt{tibercad-2.0.0-setup.exe} & Windows 32-bit     \\
    \texttt{tibercad-2.0.0-installer.bin} & Linux 32-bit self-extracting installer  \\
    \hline
    \end{tabular}
 

Linux installation procedure
------------------------------------

To install tiberCAD under Linux, download and run the self-extracting installer  ``tibercad-2.0.0_installer.bin`` 
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
and start the simulation by double clicking the input file, e.g. ``bulk.tib`` in :ref:`tut0step1` . If
not, copy the whole directory to a location in your personal area and run the examples
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
run tiberCAD by (assuming :ref:`tut0step1` )

  ``$ tibercad bulk.tib``


Bug reports / Feedback
-------------------------------------

Please send bug reports, feedback or suggestions to support@tibercad.org. When
submitting bug reports, please always include the full version number of tiberCAD you
are running. The full version number appears in the first line of output when running
the program:

::

  $ tibercad -v
  TiberCAD version 2.0.0 (x86_64-linux)

