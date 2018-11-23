.. image:: SLIMEr.png
   :align: center

|zenodo| |version|

Collection of scripts for adding SLIME (Split Lipids Into Measurable Entities) reactions into the genome-scale model of yeast. For more info, refer to the current pre-print: `SLIMEr: probing flexibility of lipid metabolism in yeast with an improved constraint-based modeling framework <https://www.biorxiv.org/content/early/2018/09/14/324863>`__.

Last update: 2018-11-22

This repository is administered by Benjamin J. Sanchez (`@BenjaSanchez <https://github.com/benjasanchez>`__), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

Installation
------------

Required Software:
~~~~~~~~~~~~~~~~~~

-  A functional Matlab installation (R2016b or higher).
-  The `COBRA toolbox for MATLAB <https://github.com/opencobra/cobratoolbox>`__.
-  For performing the lipid validation: `optGpSampler <http://cs.ru.nl/~wmegchel/optGpSampler/>`__ downloaded and unzipped in `/simulations`. Make sure the path to your solver is in your system path variable.

Dependencies - Recommended Software:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  libSBML MATLAB API (version `5.16.0 <https://sourceforge.net/projects/sbml/files/libsbml/5.16.0/stable/MATLAB%20interface/>`__ is recommended).
-  Gurobi Optimizer for MATLAB (version `7.5.2 <http://www.gurobi.com/registration/download-reg>`__ is recommended).
-  An `implementation for Matlab <https://github.com/BenjaSanchez/cmaputil/tree/master/cmaputil_matlab>`__ of the `cividis colormap <https://journals.plos.org/plosone/article/comments?id=10.1371/journal.pone.0199239>`__.

Contributors
------------

-  `Feiran Li <https://www.chalmers.se/en/staff/Pages/feiranl.aspx>`__ (`@feiranl <https://github.com/feiranl>`__), Chalmers University of Technology, Gothenburg Sweden
-  `Benjamin J. Sanchez <https://www.chalmers.se/en/staff/Pages/bensan.aspx>`__ (`@BenjaSanchez <https://github.com/benjasanchez>`__), Chalmers University of Technology, Gothenburg Sweden

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1494241.svg
   :target: https://doi.org/10.5281/zenodo.1494241
.. |version| image:: https://badge.fury.io/gh/sysbiochalmers%2Fslimer.svg
   :target: https://badge.fury.io/gh/sysbiochalmers%2Fslimer
