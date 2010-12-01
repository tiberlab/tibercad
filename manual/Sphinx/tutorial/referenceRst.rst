ReST Reference Guide

:Authors: TiberLAB srl
:Contact: support@tibercad.org

Copyright |copy| 2008, 2010 TiberCAD\ |trade|.

.. |copy| unicode:: 0xA9 
.. |trade| unicode:: U+2122

Preamble
========

This is the list of the commands that you can use, but...

.. Warning:: **NO** Tabs, only four spaces

.. note::

   Beware to the blank lines of separations between command and text.
    
   Beware to the white spaces before the inline commands and respect indentation

.. contents:: Contents


Toc (Table of Contents)
=======================

To make an automatic Index, add this code : ::

    .. contents:: Contents

This line will be replaced by a list of sections and subsections still declared in the text.

If you want to limit the depth of the list, use :  ::
    
    .. contents:: Contents
       :depth: [number]
	
*[number]* refers to the depth level of the Doctree.

To insert an anchor that binds the Section with the ToC, add :  ::

    back to Contents_
	
.. note:: 
   This is a typical linking, the anchor has the same name of the Target and ends with a **_**.

| 
| back to Contents_

Formatting
=======================

**italic**

To format a normal text in *italic*, do so : ::

	*italic*

**bold**

To format a normal text in **bold**, do so : ::

	**bold**

**code**

To format a normal text as it was ``code text``, do so : ::

	``code text``

**bulleted list**

To format a list in a bulleted way, do so : ::

    * first element
    * second element
    * third element

* first element
* second element
* third element

.. note::
   Notice the white space and the blank lines before the text
   
**numbered list**

To format a list with a numerical pointer, do so : ::

    #. first element
    #. second element
    #. third element

#. first element
#. second element
#. third element

.. note::
   Notice the white space and the blank lines before the text
   
**comment**

To insert a comment that will **not** be displayed on Output, do so : ::

    ..
      this is a comment, and to read it click on Show Source on the left

..
  Our release is one year late, this is the true

**superscript**

To format a superscrit :superscript:`text` in your documentation, do so : ::

    :superscript:`text`

**subscript**

To format a subscript :subscript:`text` in your documentation, do so : ::

    :subscript:`text`

| 
| back to Contents_


Image & Figures & Graphs
=======================

Image
-----

To display a simple image, do so : ::

    .. image:: http://www.tibercad.org/files/images/CB1r_VB1y.thumbnail.jpg
       :scale: 100 %
       :alt: This is my alt Text
       :align: center

.. image:: http://www.tibercad.org/files/images/CB1r_VB1y.thumbnail.jpg

.. note:: you can also insert a relative path like ../somewhere

| 
| back to Contents_


Figures
-------

To display a figure (image with caption), do so : ::

   .. figure:: http://www.tibercad.org/files/images/strain_trace_0.thumbnail.png
      :scale: 100 %
      :alt: Volume deformation in conic quantum dot
      :align: center
	  
      Volume deformation in conic quantum dot 

.. figure:: http://www.tibercad.org/files/images/strain_trace_0.thumbnail.png
   :scale: 100 %
   :alt: Volume deformation in conic quantum dot
   :align: center
   
   Volume deformation in conic quantum dot 
   
| 
| back to Contents_

Graphs
------

To display a simple graph of balloon type, do so : ::

   .. graph:: foo

   "start" -- "end";


.. graph:: foo

   "start" -- "end";

   
.. note::
   If you **don't see the graph**, you need the Sphinx extension based on Graphviz.
   
   To enable the extension we have to add it to the extensions list in conf.py:

| 
| back to Contents_

Math expression
===============

You then can type standard LaTeX math expressions,

**inline**
  either inline, inside text: ::

    :math:`LaTeX math expression`

:math:`a^2+b^2=c^2`


**block**
or in display mode on a new line: ::

   .. math::

      LaTeX math expressions

.. math::

   \sum_{n=0}^N x_n = y

::

   Directives (Block elements) are a general-purpose extension mechanism.
   The general syntax is:

   .. [name]:: [argument 1]
            [argument 2]
   :[option 1]: [value]

   [body]

.. note:: To enable the extension, the following line has to appear in conf.py:

   extensions = ['sphinx.ext.pngmath']
   Notice also that the sphinx.ext.pngmath extension **needs** dvipng.

**labels**

Equations are labeled with the label option and referred to using: ::

    .. math:: a^2 + b^2 = c^2
       :label: pythag

    See equation :eq:`pythag`.

.. math:: a^2 + b^2 = c^2
   :label: pythag

See equation :eq:`pythag`.

| 
| back to Contents_

Coding & Notes
==============

Coding
------

If you want to highlight a piece of code in a known format, use: ::

    .. code-block:: [language]
       :linenos:

       [body]

..
  .. code-block:: [python]
     :linenos:

     # The short X.Y version.
     version = '1.2.2'
     # The full version, including alpha/beta/rc tags.
     release = '1.2.2'

.. note:: Set **Pygment** to [language] for code highlighting.
   
   If you include **:linenos:** option, the row numbers will be displyed!

| 
| back to Contents_


Notes
-----

**note**

To create a note highlight, do so: ::

    .. note:: This is a note

.. note:: This is a note

**warning**

To create a warning highlight, do so: ::

    .. warning:: This is a warning

.. warning:: This is a warning

**see also**

To create a see also box, do so : ::

    .. seealso::

       [reST definition list]

.. seealso::
   Apples
   
   A kind of fruit.

**rubric**

Create a title **not appearing** in the table of contents by: ::

    .. rubric:: [Title]

.. rubric:: This is a Title

*Centered text*

Create a centered, boldface text block with: ::

    .. centered:: [text block]

.. centered:: Centered Text

| 
| back to Contents_

Hyperlinks & Anchors
====================

Hyperlinks
----------

**Links**

There exist two version for doing this. Either in a citation style or in an inline style :

*Citation*: ::

    This is a link_ in citation style

    .. _link: http://www.google.it

This is a link_ in citation style 

.. _link: http://www.google.it

*Inline*: ::

    In-line versions are
    `link <http://www.google.it>`_
    or `<http://www.google.it>`_
    or (autolinked) http://www.google.it.

In-line versions are
`link <http://www.google.it>`_
or `<http://www.google.it>`_
or (autolinked) http://www.google.it.


| 
| back to Contents_


Anchors
-------

**Anchor**

To define a label for any text location, precede it with: ::

    .. _[label]:

|          plus a blank line. 

.. _testA:

To reference ``[label]`` defined in any document of the project use: ::

    :ref:`[displayed text] <[label]>`

:ref:`Click here to go to TestA <testA>`

| 
| back to Contents_


Footnotes & Citations
=====================

Footnotes
---------

Footnotes:

To insert a reminder, or a note at the end of the page, create a footnote [1]_: ::

    [1]_

Once you have put the marker at the end of the page, you can insert the citation: ::

    .. [1] A footnote contains body elements, consistently indented by at least 3 spaces.
   
.. [1] A footnote contains body elements, consistently indented by at least 3 spaces.


| 
| back to Contents_


Citations
---------

Citations:

Just like a footnote, except the label is textual [CIT2002]_: ::

    In our example, we'll introduce the citation [CIT2002]_

Once you have put the marker at the end of the Document, you can insert the citation: ::

    .. [CIT2002] Just like a footnote, except the label is textual.

.. [CIT2002] Just like a footnote, except the label is textual.

| 
| back to Contents_


Markers
=======

To insert a marker that can be replaced by whatever you want, use this:  ::

    The |biohazard| symbol must be used on containers used to dispose of medical waste.

The |biohazard| symbol must be used on containers used to dispose of medical waste.

After that you can write a directive to suggest the substitution: ::

    .. |biohazard| image:: http://www.azsuperiorrestoration.com/images/biohazard.gif

.. |biohazard| image:: http://www.azsuperiorrestoration.com/images/biohazard.gif

.. note:: The syntax is :

   "|" substitution text "|" directive type "::" data directive block 

| 
| back to Contents_


Tables & Grids
=======================

Tables
------

To draw a table, depict it in this way: ::

    +--------+-------+--------+
    | Inputs |       | Output |
    +--------+-------+--------+
    |   A    |   B   | A or B |
    +========+=======+========+
    | False  |       | False  |
    +--------+-------+--------+
    | True   | False |  True  |
    +--------+-------+--------+
    | False  | True  |  True  |
    +--------+-------+--------+
    | True   |       |  True  |
    +--------+-------+--------+ 

+--------+-------+--------+
| Inputs |       | Output |
+--------+-------+--------+
|   A    |   B   | A or B |
+========+=======+========+
| False  |       | False  |
+--------+-------+--------+
| True   | False |  True  |
+--------+-------+--------+
| False  | True  |  True  |
+--------+-------+--------+
| True   |       |  True  |
+--------+-------+--------+ 

.. note:: The [**=**] divides the table Header from the content

| 
| back to Contents_


Grids
-----

You can also nest many elements and create a grid, like this: ::

    +--------+--------+-----------+
    | Header | Header with 2 cols |
    +========+========+===========+
    | A      | Lists: | **C**     |
    +--------+  * aha +-----------+
    | B::    |  * yes | | a block |
    |        |        |   of text |
    |  *hey* |  #. hi | | a break |
    +--------+--------+-----------+

+--------+--------+-----------+
| Header | Header with 2 cols |
+========+========+===========+
| A      | Lists: | **C**     |
+--------+  * aha +-----------+
| B::    |  * yes | | a block |
|        |        |   of text |
|  *hey* |  #. hi | | a break |
+--------+--------+-----------+

| 
| back to Contents_


LineBlocks
==========

To preserve breaking lines and the structure of the text, do so: ::

    | Lend us a couple of bob till Thursday.
    | I'm absolutely skint.
    | But I'm expecting a postal order and I can pay you back
    |   as soon as it comes.
    |  Ewan.

| Lend us a couple of bob till Thursday.
| I'm absolutely skint.
| But I'm expecting a postal order and I can pay you back
|   as soon as it comes.
|  Ewan.

.. note:: Notice the white space separator

| 
| back to Contents_

Bibliographic Fields
====================

To make a fields list for credits, respect this syntax: ::

    Field name :Author: author element.

**Authors**

* :Authors: authors.

**Organization**

* :Organization: organization.

**Contact**

* :Contact: contact.

**Address**

* :Address: address.

**Version**

* :Version: version.

**Status**

* :Status: status.

**Date**

* :Date: date.

**Copyright**

* :Copyright: copyright.

**Dedication**

* :Dedication: topic.

**Abstract**

* :Abstract: topic.

| 
| back to Contents_

