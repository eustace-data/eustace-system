Satellite Gridding: user guide
==============================

The satellite gridding API has been devised for processing satellite retrievals of land surface temperature: input data gets binned onto a user-defined latitute-longitude grid and then saved into a specific output directory.
Information about the various API modules can be found in  :doc:`satgrid_iris.api`.

Basic command
-------------

.. automodule:: eustace.satgrid_iris
   :members: 

Processing multiple days with Rose
----------------------------------
Two different scientific suites have been devised to process multiple sources of satellite data. 
Both suites produce satellite gridded data for a specific number of days over a specific number of years, and differ by the way multiple tasks are handled.

.. _input_parameters:

Input parameters
################

To run the analysis, a certain number of parameter values need to be assigned: this is done by opening the :code:`suite.rc` file that defines the scientific suite and changing
the values contained in the following Jinja2 templates

.. code-block:: python

  {% set OUTPATH = 'PATHTOOUTPUTDIR' %}                                          # Output data basepath
  {% set SOURCE =  'PATHTOINPUTDIR' %}                                           # Input data basepath
  {% set LST_PATTERN = 'LSTDATAFILEPATTERN' %}                                   # LST field file pattern
  {% set AUX_PATTERN = 'AUXDATAFILEPATTERN' %}                                   # AUX field file pattern

  {% set STARTYEAR = YYYYMMDD %}                                                 # Starting point: when satellite measurements have been started
  {% set ENDYEAR = YYYYMMDD %}                                                   # End point: when satellite measurements have stopped
  {% set NDAYS = 'NUMBEROFDAYS' %}                                               # Number of days to be processed in one year (default=365)
  {% set START = 'DAYSOFFSET' %}                                                 # Start day

  {% set RES = 'GRIDRESOLUTION' %}                                               # Grid resolution
  {% set NLON = 'NLOPOINTS' %}                                                   # Number of longitudinal points
  {% set NLAT = 'NLATPOINTS'' %}                                                 # Number of latitudinal points

  {% set QC_MASK_OBS = 'INT' %}                                                  # Bitmask to select the QC bits used to determine the satellite (default=1) 
  {% set QC_MASK_VALID = 'INT' %}                                                # Bitmask to select the QC bits used to determine the subset of observations 
                                                                                 # which are of interest and thought to be valid (e.g. daytime, high confidence,non-cloudy)
        									 # (default=7)

Settings for day and night analysis can be accessed through the dictionary

.. code-block:: python

  {% set SATGRID_RUN_OPTIONS = {'day': {'OUTPUT_SOURCE':'satgrid_day','DAYNIGHT':'day','QC_FILTER_OBS':0,'QC_FILTER_VALID':0},
                'night':{'OUTPUT_SOURCE':'satgrid_night','DAYNIGHT':'night','QC_FILTER_OBS':1,'QC_FILTER_VALID':1}} %}   

.. note::
   
   Satgrid_iris has predefined default values for input data directory, LST and AUX file patterns:
   
   .. code-block:: bash
   
      SOURCE = '/gws/nopw/j04/eustace/data/incoming/MODIS/'
      LST_PATTERN = '%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc'
      LST_PATTERN = '%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc'
   
   If these values corresponds to those chosen by the user, then it is not necessary to assign them. However, the script used for the analysis
   needs to be modified accordingly. With your favourite text editor, open the file 
  
   .. code-block:: bash

      /${HOME}/code/svn-eustace/research/metoffice-roses/satgrid(satgrid_rose_bunch)/app/satgrid_python/bin/launch
  
   and change the command

  .. code-block:: bash

     python2.7 -m eustace.satgrid_iris -source "${OUTPUT_SOURCE}" -qc_mask_obs=${QC_MASK_OBS} -qc_filter_obs=${QC_FILTER_OBS} -qc_mask_valid=${QC_MASK_VALID} -qc_filter_valid=${QC_FILTER_VALID} -xn ${NLON} -yn ${NLAT} -xs ${RESOLUTION} -ys ${RESOLUTION} -path ${SOURCEDIR} -pattern_lst ${LST_PATTERN} -pattern_aux ${AUX_PATTERN} -o "satgrid.${DAYNIGHT}.${TIMESTRING}.nc" "${TIMESTRING}"

  into
  
  .. code-block:: bash

     python2.7 -m eustace.satgrid_iris -source "${OUTPUT_SOURCE}" -qc_mask_obs=${QC_MASK_OBS} -qc_filter_obs=${QC_FILTER_OBS} -qc_mask_valid=${QC_MASK_VALID} -qc_filter_valid=${QC_FILTER_VALID} -xn ${NLON} -yn ${NLAT} -xs ${RESOLUTION} -ys ${RESOLUTION} -o "satgrid.${DAYNIGHT}.${TIMESTRING}.nc" "${TIMESTRING}"
      

Satgrid: processing multiple days using parametrized tasks
##########################################################

The **satgrid** suite uses `cylc <https://cylc.github.io/cylc/>`_ parametrized tasks to run several gridding procedures simultaneously,
by taking advantage of the LSF scheduler provided on CEMS. 
Here in the following are the instructions for running the suite.

1). On CEMS, in a checkout of the SVN repository, change directory to: 

.. code-block:: bash

   /${HOME}/code/svn-eustace/research/metoffice-roses/satgrid

2). Using your favourite editor, open the :code:`suite.rc` file and edit the input parameters defined in :ref:`input_parameters`.

3). Configure rose and cylc such that tasks are run on :code:`jasmin-cylc`, for more information look at :ref:`cylc-rose`.

4). Run the scientific suite by typing :code:`rose suite-run --no-gcontrol` from command line. The :code:`--no-gcontrol` option prevents cylc from launching the GUI.

5). Suite status can be monitored by using the command:

.. code-block:: bash

   ssh jasmin-cylc
   cylc monitor satgrid

Output from cylc monitor should look like

.. code-block:: bash

   satgrid_281_365 - 348 tasks                            cylc-monitor 5e12f655-7b9b-4e84-8989-0a67d8718b49
   runaheadwaitingheldqueuedreadyexpiredsubmittedsubmit-failedsubmit-retryingrunningsucceededfailedretrying
   updated: 2017-10-11T10:44:49+01                                                                         
   state summary: 43  84  213  8                                                                           
   _____________________________________________________r_u_n_n_i_n_g___t_o___s_t_o_p___a_t__20160101T0000Z
   20020101T0000Z satgrid_run_day_day281 satgrid_run_day_day282 satgrid_run_day_day283 satgrid_run_day_day284  
                  satgrid_run_day_day285 satgrid_run_day_day286 satgrid_run_day_day287 satgrid_run_day_day288 
                  satgrid_run_day_day289 satgrid_run_day_day290 satgrid_run_day_day291 satgrid_run_day_day292 
                  ........
   20030101T0000Z satgrid_run_day_day281 satgrid_run_day_day282 satgrid_run_day_day283 satgrid_run_day_day284  
                  satgrid_run_day_day285 satgrid_run_day_day286 satgrid_run_day_day287 satgrid_run_day_day288 
                  satgrid_run_day_day289 satgrid_run_day_day290 satgrid_run_day_day291 satgrid_run_day_day292 
                  ........

**satgrid** will store output data into the directory defined in :ref:`input_parameters`. It also will check for missing output files, storing their corresponding labels into a :code:`data_YYYY.txt` file.

.. warning:: 

   it is recognised that Rose chokes when handling with a large number of tasks per cycle (e.g. :code:`NTASKS>200`). A large number of tasks could affect the polling procedure that Rose adopts to check the status of running jobs. 
   This in turn slows down the submission (or the interruption) of waiting (running) tasks [2]_. If a large number of days per year needs to be processed, then it is better to adopt an incremental strategy (e.g. by dividing the number of tasks per year into smaller bunches and executing a specific suite definition for each bunch).

Satgrid_rose_bunch: processing multiple days with Rose Bunch
############################################################

The **satgrid_rose_bunch** suite uses the built-in application `Rose Bunch <http://www-nwp/~fcm/rose-doc/rose/rose-rug-task-run.html#rose-task-run.util.rose_bunch>`_ 
for running multiple command variants in parallel under a single job. 
Using this suite allows to run a large number of tasks per cycle, without being bothered by the issues described in the previous section. It also reduces the number of 
output information printed by :code:`cylc monitor`, allowing users to easily check tasks status. 
However, since a limited amount of tasks can be assigned to a single job, it usually takes more time to process all the data.
Users are advised to use this suite when data over a large range of years need to be processed.

Here in the following are the instructions for running the suite.

1). On CEMS, in a checkout of the SVN repository, change directory to: 

.. code-block:: bash

   /${HOME}/code/svn-eustace/research/metoffice-roses/satgrid_rose_bunch

2). Using your favourite editor, open the :code:`suite.rc` file and edit the input parameters defined in :ref:`input_parameters`.

.. warning:: 

   under load balancing systems such as PBS, Slurm or LSF, it is necessary to set resource requests to reflect the resources required by running multiple commands at once.
   E.g. if one command would require 1GB memory and you have configured your app to run up to 4 commands at once then you will need to request 4GB of memory.
   Moreover it is important to modify the wall time and CPU computation time values, to ensure that all the data are processed before the job will espire.
   This can is done by properly changing the values of the following parameters
   
   .. code-block:: python
   
      [[[job submission]]]                
      method = SCHEDULERNAME            # scheduler
      execution time limit = PTXXX      # wall-time, e.g. XXX=30H = 30 hours, XXX=10M = 10 minutes 
              [[[directives]]]
	      -q = QUEUENAME            # queue system
	      -c = HH:MM                # CPU-time
              -R = "rusage[mem=XXXX]"   # required total memory per job, e.g. XXXX=1000 corresponds to 1GB of memory

3). Configure rose and cylc such that tasks are run on :code:`jasmin-cylc`, for more information look at :ref:`cylc-rose`.

4). Run the scientific suite by typing :code:`rose suite-run --no-gcontrol` from command line. The :code:`--no-gcontrol` option prevents cylc from launching the GUI.

5). Suite status can be monitored by using the command:

.. code-block:: bash

   ssh jasmin-cylc
   cylc monitor satgrid

Output from cylc monitor should look like

.. code-block:: bash

   satgrid_rose_bunch_130_180 - 31 tasks                  cylc-monitor 34cd5b1a-4e62-4522-bc14-c2445a9695f7
   runaheadwaitingheldqueuedreadyexpiredsubmittedsubmit-failedsubmit-retryingrunningsucceededfailedretrying
   updated: 2017-10-11T10:00:17+01                                                                         
   state summary: 28  3                                                                                    
   _____________________________________________________r_u_n_n_i_n_g___t_o___s_t_o_p___a_t__20160101T0000Z
   20030101T0000Z satgrid_run_day satgrid_run_night                                                        
   20040101T0000Z satgrid_run_day satgrid_run_night                                                        
   20050101T0000Z satgrid_run_day satgrid_run_night                                                        
   20060101T0000Z satgrid_run_day satgrid_run_night                                                        
   20070101T0000Z satgrid_run_day satgrid_run_night                                                        
   20080101T0000Z satgrid_run_day satgrid_run_night                                                        
   20090101T0000Z satgrid_run_day satgrid_run_night                                                        


**satgrid_rose_bunch** will store output data into the directory defined in :ref:`input_parameters`. It also will check for missing output files, storing their corresponding labels into a :code:`data_YYYY.txt` file.

.. _cylc-rose:

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

This will tell Rose to run tasks only on :code:`jasmin-cylc` host [1]_.

.. [1] Due to task-job communication issues experienced on :code:`cems-sci1.cems.rl.ac.uk`, it has been decided to rely only on the :code:`jasmin-cylc.ceda.ac.uk` host.
.. [2] Moreover, checking tasks status by using :code:`cylc monitor` becomes quite unpractical.
