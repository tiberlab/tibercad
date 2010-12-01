
**Tutorial Name**


:Authors: TiberLAB srl
:Contact: support@tibercad.org

Copyright |copy| 2008, 2010 TiberCAD\ |trade|.

.. |copy| unicode:: 0xA9 
.. |trade| unicode:: U+2122

.. contents:: Contents

Preamble
============

This document provides an overview of the features of a Multiscale Simulator.

Regarding the language, almost all of the Reference Guide keywords are supported (those marked with an asterisk are
not fully supported yet):

::
 
  module, import, extends, template, query, public, private, protected, guard, init,
  overrides, each, before, after, for, if, elseif, else, let, elselet, trace*, macro*,
  file, mode, text_explicit*, code_explicit*, super, stdout

For all the details about these keywords and the TiberCAD multiscale simulator, you can consult
`the official website <http://www.tibercad.org/>`_ available for the TiberCAD's project.

| 
| back to Contents_

Development
============

The TiberCAD's development team does its best to make periodical bug free release. 

    Following is the roadmap chart :

+---------+---------------------------------+---------------------------------+---------------------------------+
| Early   |         TiberCAD 1.2.2          |    Next Release TiberCAD 1.3    | Commercial Release TiberCAD 2.0 | 
+=========+=================================+=================================+=================================+
| Free    | OK                              | OK                              |                                 |
+---------+---------------------------------+---------------------------------+---------------------------------+
| Fee     |                                 |                                 | OK                              |
+---------+---------------------------------+---------------------------------+---------------------------------+

The next release will come at March 2011.

|

    The main commercial release is scheduled for August 2011 [#]_.

| 
| back to Contents_


Modules and Projects
============================

TiberCAD is a multiscale CAD tool for the simulation of modern nanoelectronic and optoelectronic devices. 
Possible applications of TiberCAD range from nanoelectronics  to laser technologies including molecular electronics and bio-devices.


    .. image:: http://www.tibercad.org/files/images/CB1r_VB1y.thumbnail.jpg
	
.. figure:: http://www.tibercad.org/files/images/strain_trace_0.thumbnail.png
   :scale: 100 %
   :alt: Volume deformation in conic quantum dot
   :align: center

   Volume deformation in conic quantum dot 
   [leave a blank line]

| 
| back to Contents_

Physical models 
===============

Some of the features are detailed below [#]_:

* Strain/stress modelization, including pyro- and piezoelectric effects, non-linear strain, converse piezoelectric, external forces
* Classical (Drift-Diffusion, hydrodynamic) particle transport and  Poisson
* calculation, Quantum Current calculation (quantum drift-diffusion,
* Non-Equilibrium Green Function)
* Electrons, holes and excitons dynamic
* Heat balance model: electron and  hole  dissipation, microscopical  heat  model
* Quantum physics for  continuous media, including Envelope  Function Approximation, k.p theory.
* Atomistic quantum description, including Empirical Tight-Binding model and ab-initio methods
* Electromagnetic solver for optical fields

  .. image:: http://www.tibercad.org/files/images/mesh_dot.thumbnail.png

| 
| back to Contents_

Math expressions
---------------------------------------------

You can insert the Pythagoras formula :  :math:`a^2+b^2=c^2`

|
|    or make a complex expression with :
	
.. math:: 

   \sum_{n=0}^N x_n = y

| 
| back to Contents_

Get TiberCAD
------------

A  free trial version of TiberCAD  is now available for  `download <http://www.tibercad.org/features_available_version_122>`_.


| 
| back to Contents_

Updates
=======
.. note::

   Please note  that   not  all  the  features  of  the full  Tibercad  package  suite  are  available  in this  trial  version.
   The  first  commercial  release of  Tibercad will be  issued  soon. (see http://www.tiberlab.com)

| 
| back to Contents_

.. [#] The deadline may vary
.. [#] Further features can be added