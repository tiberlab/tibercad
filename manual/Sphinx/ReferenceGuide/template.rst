..       replace [ModuleName]

..   <marker>

.. index:: double:Getting Started;[ModuleName]
.. _[ModuleName]Getting:

ModuleName (Main Title)
=================================================

..  Use Bold (**text**) and Italic (*text*) for emphasis

Overview (Paragraph Level)
----------------------------

[Lorem **ipsum dolor** sit amet, *consectetur* adipiscing elit. Duis et risus a magna auctor accumsan. 
Aenean ac lobortis nulla. Vestibulum eget turpis id ante varius laoreet. 
Sed mi turpis, tristique sit amet sagittis vitae, laoreet placerat justo. Maecenas lacinia mi arcu.] 

Description
----------------------

..  Use Backquotes ``keyword`` for keywords and filenames

[Lorem ``ipsum dolor`` sit amet, consectetur adipiscing elit.]

Mini theoretical intro (Sub-Paragraph)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..  to insert math inline use it :math:`latex inline` .

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. :math:`\phi [cm^-1]`
Duis et risus a magna auctor accumsan. Aenean ac lobortis nulla. Vestibulum eget turpis id ante varius laoreet. 
Sed mi turpis, tristique sit amet sagittis vitae, laoreet placerat justo.]

.. math::
   :nowrap:
   :label:
   
   \[
   \frac{\partial}{\partial x_j}C_{ijlk}\frac{\partial u_l}{\partial x_k} = f_i
   \]

Website linking 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..  HTML links are generated automatically, but if you want to change test, use : `Link text`_ and at the end of file insert the reference

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. 
Fusce id massa sit amet nibh pellentesque aliquam.] `Fastlink`_

Note Block
^^^^^^^^^^^^^^^^^^^^^^^^^^

..  this inserts a Note Block, image can be changed. See the end of File.

 
|warn| 
       [Lorem ipsum dolor sit amet, consectetur adipiscing elit.
       Pellentesque vulputate ultrices odio, quis dictum purus commodo eu.
       Fusce id massa sit amet nibh pellentesque aliquam.]


Footnotes
^^^^^^^^^^^^^^^^^^^^^^^^^^

..  to insert a Footnote, use this syntax [#]_ and at the end of file insert the reference in the rubric

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. [#]_ 
Fusce id massa sit amet nibh pellentesque aliquam. [#]_]

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. 
Fusce id massa sit amet nibh pellentesque aliquam.]

Figures
^^^^^^^^^^^^^^^^^^^^^

..  to insert a figure with caption use the follow example

.. figure:: ../data/elasticity01.png
   :align: center
   :scale: 50%

   Caption (with a blank line)


..  to insert an anchor of a figure use a representive name preceeded by "_" and followed by ":" with no text

.. _fig.elasticity01: 

..  to insert a line break use ---- (four minus)

----
  
[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. Fusce id massa sit amet nibh pellentesque aliquam. Nunc vel est eget urna sodales scelerisque. Vivamus blandit euismod lacinia. Aliquam dolor mi, semper a porttitor malesuada, eleifend nec ligula]

Click to return to the :ref:`first figure of elasticity<fig.elasticity01>`


Example
-------

..  Use the :: operator to highlight Code-Blocks

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. Fusce id massa sit amet nibh pellentesque aliquam. Nunc vel est eget urna sodales scelerisque. Vivamus blandit euismod lacinia. Aliquam dolor mi, semper a porttitor malesuada, eleifend nec ligula]

::

  Device
    {
     dimension = 3 
     meshfile = mesh.msh
     mesh_units = 1e-9
     material = AlGaN
     x = 0.2
     x-growth-direction = (-1,0,1,0) 
     y-growth-direction = (-1,2,-1,0) 
     z-growth-direction = (0,0,0,1)
     Region Well  {material = GaN}
    }


Keywords
-------------------


Tables
^^^^^^^^^^^^^^^^^^^

..  instead of using the table depiction method, insert your Latex table as Math Block

[Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque vulputate ultrices odio, quis dictum purus commodo eu. Fusce id massa sit amet nibh pellentesque aliquam. Nunc vel est eget urna sodales scelerisque. Vivamus blandit euismod lacinia. Aliquam dolor mi, semper a porttitor malesuada, eleifend nec ligula]

..  math::
    :nowrap:
    :label:

    \begin{table}[!h]
    \center
    \begin{tabular}{l|c|l|l|r}
    \multicolumn{5}{c}{\textbf{Keywords Table}} \\
    \hline
    \textbf{Keyword} & \textbf{Type} & \textbf{Description} & \textbf{Unit} & \textbf{Default} \\
    \hline
    \hline
     &  &  &  &  \\
    \texttt{keyword01} & Type01 & Desc01 & $[cm^{-1}]$ & Default01 \\
    \texttt{keyword02} & Type02 & Desc02 & $[cm^{-1}]$ & Default02 \\
    \texttt{keyword03} & Type03 & Desc03 & $[cm^{-1}]$ & Default03 \\
    \texttt{keyword04} & Type04 & Desc04 & $[cm^{-1}]$ & Default04 \\
    \hline
    \end{tabular}
    \caption{Keyword listing for [ModuleName]}
    \end{table}

Keyword Explanation
^^^^^^^^^^^^^^^^^^^^^

| ``keyword01``
|       Explanation of keyword01
| 
|       Further explanations

| ``keyword02``
|       Explanation of keyword02
| 
|       Further explanations



Run TiberCAD
-----------------------

..  Using | you can force the layout of text similar at Source Code

Now we can run TiberCAD
      | by double clicking on bulk.tib file (in Windows)
      | 
      | or by command line in linux: tibercad bulk.tib


..   </marker>

..  _Fastlink: http://www.tibercad.org

..   rubric:: Footnotes

..   [#] text of first Footnote
..   [#] text of second Footnote

.. |more| image:: ../data/more.png
    :scale: 50%

.. |warn| image:: ../data/warn.png
    :scale: 50%

.. |idea| image:: ../data/idea.png
    :scale: 50%
