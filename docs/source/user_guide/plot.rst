Plot Tab
========

The **Plot** tab provides a built-in plotting utility for visualizing cross
sections, S-factors, and calculation results. This tab requires the QWT library
to be enabled at compile time (``-DUSE_QWT=ON``).

.. note::

   The built-in plotter is intended as a quick visualization tool during
   analysis. For publication-quality figures, use the output files
   (``AZUREOut_aa=*_R=*.out`` and ``AZUREOut_aa=*_R=*.extrap``) with an
   external plotting program.

Selecting Segments
------------------

The tab displays two lists of segments, corresponding to the segments defined
in the **Segments** tab:

- **Segments From Data** -- data-driven calculation segments.
- **Segments Without Data** -- pure calculation segments.

Click a segment to select it (the background color changes). Click again to
deselect. Multiple segments can be selected simultaneously.

Drawing Plots
-------------

After selecting one or more segments, click the **Draw** button in the
lower-left corner to render the plot.

- **R-matrix calculations** are shown as red lines.
- **Data points** are shown in black with different marker shapes for different
  segments.

Axis Configuration
------------------

X-Axis
^^^^^^

Choose the horizontal axis quantity:

- **CoM Energy** -- center-of-mass energy
- **Excitation Energy** -- excitation energy of the compound nucleus
- **CoM Angle** -- center-of-mass angle

Y-Axis
^^^^^^

Choose the vertical axis quantity:

- **Cross Section** -- cross section (barns or barns/sr)
- **S-Factor** -- astrophysical S-factor (MeV b or MeV b/sr)

Both axes have a **Log** checkbox to switch to logarithmic scale.

The axis title follows what the selected segments actually plot, rather than
the radio button alone. A differential segment is labelled ``[b/sr]`` and an
angle-integrated one ``[b]``; a phase-shift segment is labelled in degrees; and
the two polarization observables are named as themselves --
**Analyzing Power**, which is dimensionless, and
:math:`P\,d\sigma/d\Omega` in ``[b/sr]``. If a selection mixes quantities
whose units differ, the title carries no units rather than asserting the wrong
ones.

.. note::

   Drawing an analyzing-power or **Polarization x Cross Section** segment
   switches the y-axis to a linear scale automatically, overriding both
   settings. This is a necessity rather than a preference: both observables
   take negative values -- :math:`A_y` is a dimensionless ratio negative over
   much of its range, and :math:`P\,d\sigma/d\Omega` is negative wherever the
   polarization is -- so a logarithmic axis cannot represent either, and an
   S-factor conversion has no meaning for them. Negative points are drawn; for
   every other observable the non-positive values that a logarithmic axis
   cannot show are still filtered out.

.. tip::

   Plotting with **Excitation Energy** on the x-axis is very useful for
   comparing data from different entrance and exit particle channels, since
   all channels share the same excitation energy scale.

Interacting with Plots
-----------------------

**Zoom in**
   Click and hold the left mouse button, then drag to highlight a region.

**Zoom out**
   Right-click anywhere on the plot.

**Pan**
   Hold the center mouse button and drag to shift the view.

When first plotted, axes are automatically scaled to fit all data.

Exporting and Printing
----------------------

**Export...**
   Save the current plot to an image file. Most standard image formats are
   supported.

**Print...**
   Send the current plot to a printer.

Both export and print use the currently displayed view, so you can zoom in
before exporting to customize the plotted area.
