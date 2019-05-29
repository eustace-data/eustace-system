Satstace guide
===============================

This guide shows how to run the Satstace system for the production of air temperature estimates (with estimates of uncertainty) for all surfaces of earth.
Very few simple instructions have to be followed, as explained below:

1). On CEMS, in a checkout of the SVN repository, change directory to: 

.. code-block:: bash

   /${HOME}/code/svn-eustace/system/src/roses/satstace

2). If you are running tests, then make a new subdirectory on CEMS in 

.. code-block:: bash

   /work/scratch/${USERNAME}

where :code:`${USERNAME}` has to be replaced by your own username [1]_ . If you are running the production of output, then you should create a folder in

.. code-block:: bash

   /gws/nopw/j04/eustace/data/internal/D2.2/${REVISIONNUMBER}/${DATE}

where :code:`${REVISIONNUMBER}` refers to the version of the system adopted, and :code:`${DATE}` refers to the date the new data set has been created. :code:`${REVISIONNUMBER}` should
have the following standard format

.. code-block:: bash

   ${REVISIONNUMBER}=RXXXX

where :code:`XXXX` is the system version number.

3). Inside :code:`/${HOME}/code/svn-eustace/system/src/roses/satstace`,have a look at the available json descriptors

.. code-block:: bash

   $ls *.json
   ice.json land.json ocean.json

These descriptors encode information about the input sources needed by SATSTACE to run the analysis and their corresponding path: before running the analysis, you should check all the information
contained in the descripors is correct.

4). Using your favourite editor, open the :code:`environment.sh` bash script and perform the following changes
    * set :code:`OUTPUTDIRECTORY` to the directory name defined at step 2)
    * :code:`export VERSIONCONTROL=--allow_unversioned_code`
    * :code:`export CHECKSUM=--nochecksum`
    
5). Execute the :code:`makecatalogue_inputs.sh` bash script: this will call :doc:`eumopps` to build the catalogue of relevant inputs and operations needed to run SATSTACE.

.. code-block:: bash

   ./makecatalogue_inputs.sh

If the script runs properly, it should produce a `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ catalogue inside :code:`OUTPUTDIRECTORY`: if you want to inspect catalogue content
just use the tools provided by :doc:`eumopps`.

6). Using your favourite editor, open the :code:`suite.rc` file and modify the value of the :code:`CATALOGUE` variable
 
.. code-block:: python

   {% set CATALOGUE = '/work/scratch/${USERNAME}/catalogue.nc' %}

by setting :code:`${USERNAME}` to your own username.

7). Configure rose and cylc such that tasks are run on :code:`jasmin-cylc`, for more information look at :ref:`cylc-rose1`.

8). Run the scientific suite by typing :code:`rose suite-run --no-gcontrol` from command line. The :code:`--no-gcontrol` option prevents cylc from launching the GUI.

9). Suite status can be monitored by using the command:

.. code-block:: bash

   ssh jasmin-cylc
   cylc monitor satstace

Satstace output will be stored into the working directory defined at step 2).

Tips and tricks
---------------

Sometimes the cluster computing facilities we use do not work seamlessly: sometimes jobs could fail, we could run out scientific software licenses, or adopted analysis frameworks do not run properly.
In this case, issues could be mitigated by taking proper actions: here in the following you will find how to mitigate some of the issues encountered when running SATSTACE using `JASMIN <http://www.ceda.ac.uk/projects/jasmin/>`_ facilities.

1). **Running out of** `IDL <https://www.harrisgeospatial.com/SoftwareTechnology/IDL.aspx>`_ **licenses**: at the core level, regression analysis for land and ocean surfaces is performed by 
`IDL <https://www.harrisgeospatial.com/SoftwareTechnology/IDL.aspx>`_ procedures, which are then wrapped by python code to be integrated into the system. Unfortunately `IDL <https://www.harrisgeospatial.com/SoftwareTechnology/IDL.aspx>`_
is not free, and a limited amount of its instances can be run, depending on the number of available licenses.
When too many instances run together, some of their corresponding jobs will fail, due to the lack of available licenses: the rose suite :code:`suite.rc` has been edited to take into account of this issue.
However, further mitigation can be obtained by modifying the following field
 
.. code-block:: python

   {% set CHUNKSIZE = 5 %}

which sets the maximum number of concurring taks in a given cycle.

2). **I/O issues when too many jobs are running**: `cylc <https://cylc.github.io/cylc/>`_ does not handle scientific suites that run many parallel tasks at the same time. Moreover, issues have been encountered
when running jobs on `JASMIN <http://www.ceda.ac.uk/projects/jasmin/>`_ facilities, resulting in having some of the jobs disappeared along with their expected output. This also caused `cylc <https://cylc.github.io/cylc/>`_ not to be able
to complete its cycles, as it did not recognized that some of the submitted tasks have disappeared. This phenomenon becomes more evident when one tries to run the regression analysis for all the surfaces at the same time.
To mitigate this issue, the following approach should be followed

  a). Run the analysis for each single surface at time: this can be done by editing the rose suite and commenting all the tasks but one 
  
  .. code-block:: python

	      [[[R/P1/{{(NUM_SATSTACE_LAND/CHUNKSIZE)|roundup}}]]]
		  graph = """satstace_land_subtasks"""

      #        [[[R/P1/{{(NUM_SATSTACE_ICE/CHUNKSIZE)|roundup}}]]]
      #           graph = """satstace_ice_subtasks"""

      #        [[[R/P1/{{(NUM_SATSTACE_OCEAN_AATSR/CHUNKSIZE)|roundup}}]]]
      #           graph = """satstace_ocean_AATSR_subtasks"""

  b). Tune the number of simultaneously active cycles
    
  .. code-block:: python

    {% set MAX_SIMULTANEOUS = 30 %}

  c). Use an incremental strategy: instead of processing the entire block of input data for a given surface, divide it in smaller blocks and run the analysis in sequence, one block at time.
  This can be done by tuning the following fields
    
  .. code-block:: python

     initial cycle point = 15
     final cycle point = 115

  inside the rose suite. It is also useful to identify rose suites with the name of the task performed, the initial and final cycle points values: e.g. given the following setup
 
  .. code-block:: python

	      [[[R/P1/{{(NUM_SATSTACE_LAND/CHUNKSIZE)|roundup}}]]]
		  graph = """satstace_land_subtasks"""

      #        [[[R/P1/{{(NUM_SATSTACE_ICE/CHUNKSIZE)|roundup}}]]]
      #           graph = """satstace_ice_subtasks"""

      #        [[[R/P1/{{(NUM_SATSTACE_OCEAN_AATSR/CHUNKSIZE)|roundup}}]]]
      #           graph = """satstace_ocean_AATSR_subtasks"""
      ...
    initial cycle point = 15
    final cycle point = 115
      
  one could submit the rose suite with the :code:`--name` flag and name it as :code:`satstace_land_15_115`

  .. code-block:: bash

     rose suite-run --no-gcontrol --name=satstace_land_15_115
  
  This would make suite traceability easy.
  
  d.) Periodically check suite status and the amount of produced output inside :code:`OUTPUTDIRECTORY`.

.. _cylc-rose1:

Configuring cylc and Rose
-------------------------

Using your favourite editor, open the :code:`~/.cylc/global.rc` file and edit the following lines

.. code-block:: bash

   [editors]
   terminal = emacs
   gui = emacs

   [suite host scanning]
   hosts = jasmin-cylc

   [communication]
   method=http

This will force Rose to scan taks only on :code:`jasmin-cylc` host, and to use the http method for communicationg with them.
Then open :code:`~/.metomi/rose.conf` and edit the following lines

.. code-block:: bash

   [rose-suite-run]
   hosts=jasmin-cylc

This will tell Rose to run tasks only on :code:`jasmin-cylc` host [2]_.


.. [1] More specific paths can be adopted, depending on users' choice, e.g. :code:`${USERNAME}/satstace/`, or :code:`${USERNAME}/${NAME1}/.../${NAMEN}`
.. [2] Due to task-job communication issues experienced on :code:`cems-sci1.cems.rl.ac.uk`, it has been decided to rely only on the :code:`jasmin-cylc.ceda.ac.uk` host.