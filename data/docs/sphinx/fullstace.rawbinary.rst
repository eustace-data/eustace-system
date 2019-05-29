Full Raw binary generation guide
================================

This guide describes how to run the Fullstace system to produce rawbinary output.

Instructions:


1). On CEMS, in a checkout of the SVN repository, change directory to: 

.. code-block:: bash

   /${HOME}/code/svn-eustace/system/src/roses/fullstace/rawbinary


2).Confirm if the current analysis run is for testing or production purposes 


3). View the json descriptor files to confirm correct settings for input data, date ranges and quality control flag settings

.. code-block:: bash

   $ls *.json
   rawbinary_inputs.json rawbinary_generation.json 

These descriptors encode information about the input sources needed by FULLSTACE to run the analysis and their corresponding path: before running the analysis, you should check all the information
contained in the descriptors is correct.


4). To run tests for only insitu_land input data, with optional qc flags settings, follows these steps

Change directory to the test subdirectory

.. code-block:: bash

   cd test

Execute one or more of the following scripts to create rawbinary files for only insitu land inputs with or without qc flags used to reduce the number of valid input sources


.. code-block:: bash

   source test_create_rawbinary_files_land_noqc_ho.sh
   source test_create_rawbinary_files_land_qc_ho.sh

These tests use the .json files within this test subdirectory to produce output quickly for a single year.  


5). To create a full production run of the rawbinary files, follow these steps.

Change directory to 

.. code-block:: bash

   /${HOME}/code/svn-eustace/system/src/roses/fullstace/rawbinary

Execute the :code:`create_rawbinary_files.sh` bash script: this will call :doc:`eumopps` to build the catalogue of relevant inputs and operations needed to run FULLSTACE.

.. code-block:: bash

   source create_rawbinary_files.sh

This script can be modified directly, or its component commands can be run in sequence if alternative output directories, or other variations in settings are required. 


If the script runs correctly, it should produce a `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ catalogue inside the users directory at :code:`/work/scratch/$USER/rawbinary`: if you want to inspect catalogue content use the tools provided by :doc:`eumopps`.

Once verified and considered of production quality, copy the ouput to an appropriate directory based on recent code version and the analysis run date (replace revision and date). ' 
e.g. 

:code: `gws/nopw/j04/eustace/data/internal/D2.2/$REVISION/$DATE`:


6). Troubleshooting 

In production use, with all data source inputs, the above script will launch many jobs on the LSF scheduler.

Progress of the jobs can be checked with the commands

.. code-block:: bash

   bhist
   bjobs

Text files with any error or output information that had been sent to stderr or stdout are produced in 

:code: `/work/scratch/$USER/lsf`

If there are problems or you launch the script in error, you can  use 

.. code-block:: bash

   bkill 0

to kill every LSF job you have running including those ones (use with caution!)
