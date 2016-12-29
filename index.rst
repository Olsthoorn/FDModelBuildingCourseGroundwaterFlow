.. mfpy documentation master file, created by
   sphinx-quickstart on Tue Oct 18 23:05:44 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Finite Difference Grounwater Modeling in Python
============================================================

.. image:: ./_figures/ThreeWellPlusInactive.png
  :align: center
  :scale: 50%

::

    Finite Difference Models are derived and implemented completely in Python.
    The theory and construction of these models can be used in their own right
    or may serve as a thorough introduction in groundwater modeling with available
    codes especially with MODFLOW, MT3DMS, MODPATH and SEAWAT.
    At the end of this course we have built from ground on a powerful 3D steady state and transient finite difference groundwater code completely in Python functions and also a powerful 3D particle tracking funcion capable of tracking millions of particles simultaneously. We also have seen the versatile use of this code. The finite difference model functions are compatible with MODFLOW and MODPATH. The only limitation is that the finite difference functions allow just fixed-head and prescribed flow boundaries. This is to limit clutter and keep the finite-difference model functions in this course as lean and simple to use as possible, but the full theory is presented in the first chapter. This limitation is not a problem in most cases as, it's easy to model so-called general-head boundaries by setting appropriate conductivities.
    For more advanced finite difference modeling and use of more specific packages, one should use the USGS codes MODFLOW, MODPATH, MT3DMS and SEAWAT directly. A Matlab interface as developed and used by me and my students for the last 8 years in many projects is available under project mfLab on SourceForge.org. A Python interface is available under the name PyFlow as developed by the USGS.
.. sectionauthor:: Prof. dr.ir. Theo Olsthoorn (emer. TU Delft, the Netherlands)


.. toctree::
   :maxdepth: 4
   :caption: Contents:

01_Numerical_grw_modeling.ipynb
02_fin_dif_modeling.ipynb
03_fdm_as_python_func.ipynb
04_computing_flows.ipynb
05_grid_class.ipynb
06_axial_symmetry.ipynb
07_stream_lines.ipynb
08_transient_flow.ipynb
09_particle_tracking.ipynb


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
