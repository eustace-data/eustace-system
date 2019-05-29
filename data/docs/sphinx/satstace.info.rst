Satstace: how does it work?
===============================

Satstace is a sub-part of the EUSTACE end-to-end system.
It provides surface air temperature estimates (with estimates of uncertainty) for all surfaces of earth, derived from satellite surface skin temperature retrievals.
The whole production procedure is characterized by three different steps:

1). *Surface skin temperature retrievals*. 
 
Values of surface skin temperature are collected from satellites and insitu sources, for each day since 2000. On CEMS, data can be found at:

.. code-block:: bash

   /gws/nopw/j04/eustace/data/incoming/
   /gws/nopw/j04/eustace/data/internal/

Data are usually stored into NetCDF files.

2). *Data preprocessing*. 

Satellite data are aggregated onto a regular latitude-longitude grid. For ice and ocean temperature retrievals, this procedure is already performed at step 1), while land surface data are preprocessed by using the :doc:`satgrid` utility.

3). *Surface air temperature fields production*. 

For each surface dataset, a specific regression algorithm is applied to determine the functional relation between skin and surface air temperature. The algorithm is trained by using skin temperature retrievals and insitu measurements of surface air temperature as the corresponding input and output variables.
On a checkout of the SVN repository, the implementation of the regression algorithms can be inspected at
 
.. code-block:: bash

   /${HOME}/code/svn-eustace/system/src/eustace/surfaceairmodel

The regression analysis can be (potentially) performed simultaneously for all the different surfaces, by running few bash scripts and a `Rose <http://www-nwp/~fcm/rose-doc/rose/rose.html>`_ suite, as explained in :doc:`satstace.guide`.
