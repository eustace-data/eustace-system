Analysis System Guide
=====================

This page shows users how to run the eustace end-to-end system to perform the advanced standard method analysis.

Json descriptors
----------------

All the relevant information about the inputs to be processed, operations that need to be performed and the output that has to be produced
are encoded into json descriptors. 
The latter are processed by :ref:`EUMOPPS <eumopps>` to generate a detailed catalogue used to run the analsysis.

In this case all the relevant information are stored into two different json files: if you have checked out a working copy of
the system into your :code:`${HOME}` directory, they can be found in

.. code-block:: console

   $ls ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/*json
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/advstd.json
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/advstd_inputs.json

If you look in :code:`advstd_inputs.json`, this is what you should see

.. code-block:: json

    {
	"python_class" : "eumopps.catalogue.catalogue.Catalogue",

	"datasets" : [

	{
	    "python_class" : "eumopps.catalogue.dataset.CatalogueDataSet",
	    "name"  : "locationlookup_insitu_land",
	    "path"	: "/work/scratch/eustace/rawbinary",
	    "subsets" : [
	    {
		"python_class" : "eumopps.catalogue.dataset.CatalogueDataSubset",
		"layout" : 
		{
		    "python_class" : "eumopps.catalogue.storage.DataStorageFiles",
		    "patterns" : [ "locationlookup_insitu_land.bin" ]
		}	    
	    }
	    ]
	},

	{
	    "python_class" : "eumopps.catalogue.dataset.CatalogueDataSet",
	    "name"  : "locationlookup_satellite",
	    "path"	: "/work/scratch/eustace/rawbinary",
	    "subsets" : [
	    {
		"python_class" : "eumopps.catalogue.dataset.CatalogueDataSubset",
		"layout" : 
		{
		    "python_class" : "eumopps.catalogue.storage.DataStorageFiles",
		    "patterns" : [ "locationlookup_satellite.bin" ]
		}	    
	    }
	    ]
          },

This descriptor contains information about all the relevant inputs the system has to process, like in situ data sources, satellite data produced by :ref:`SATSTACE <satstace_guide>`, correlation ranges etc.
These data sources come into the form of rawbinary files, and their generation should precede the infilled analysis.

On CEMS, these intermediate data sources should be stored in a directory like :code:`/work/scratch/eustace/rawbinary`.
If such directory exists, then rawbinary files have already been produced, otherwise you should have a look at :ref:`raw_binary_files_production`.

Inside the other json descriptor, :code:`advstd.json`, you should recognize different key sections like

.. code-block:: json

   {
      "python_class" : "eumopps.catalogue.catalogue.Catalogue",

      "operations" : [

      {
	  "python_class" : "eumopps.catalogue.operation.Operation",
	  "name" : "climatology_input",
	  "runmodule" : 
	  {

and

.. code-block:: json

   {"step" : { "python_class" : "eumopps.catalogue.step.StepDaily", "start" : "20060101000000", "end" : "20061231000000" },
        "newdatasets": [
        {
	  "python_class" : "eumopps.catalogue.dataset.CatalogueDataSet",
	  "name"  : "measurement_climatology",
	  "subsets" : [

and also

.. code-block:: json

    {
	"python_class" : "eumopps.catalogue.operation.Operation",
	"name" : "climatology_solve",
	"runmodule" : 
	{
	    "python_function" : "eustace.analysis.advanced_standard.examples.example_eustace.solve",
	    "storage_climatology" :
	    {

The first section describes a specific operation, the one that will read data from rawbinary files sources and build the measurement input used for the infilled analysis.
If you look in detail into the descriptor, you will notice sections describing similar operations, but for a different model component, like

* :code:`measurement_climatology`
* :code:`measurement_large_scale`
* :code:`local_input_and_solve`

The names assigned to each operation are arbitrary, and can be changed by users if necessary. What should not be changed are the other input information of this section, like the name of the python module called by :ref:`EUMOPPS <eumopps>` and its corresponding input.

The second section represents the time interval the system should cover when performing the analysis: the start and end date can be changed by users to analyze different time periods. It can be found
where the production of inputs for climatology, large and local scale models is specified.

The last section represents operations that solve each model component of the system: an exception is done for the local model, where the computation of inputs and the model solution are performed together.

At the bottom of the json file you will find the section labeled as :code:`output_grid`: it describes the final output production onto a rectilinear grid.
Also in this case users can specify the time period they want to cover when producing the output.

Tuning the climatology model
----------------------------

Inside :code:`advstd.json`, you will find the following field

.. code-block:: json

    {"covariates_descriptor" : "/gws/nopw/j04/eustace/data/internal/climatology_covariates/covariates.json",
    }

occurring in several places. This field refers to another json descriptor, used for describing the shape of the climatology model.

The json descriptor should look like

.. code-block:: json

    {"latitude_harmonics":{"element":{"python_class":"eustace.analysis.advanced_standard.elements.latitudeharmonics.LatitudeHarmonicsElement"},
			  "hyperparameters":{"hyperparameter_1":-1.15,
					      "hyperparameter_2":-1.15,
					      "hyperparameter_3":-1.15,
					      "hyperparameter_4":-1.15}},
    "altitude":{"element":{"python_class":"eustace.analysis.advanced_standard.elements.geography_based.GeographyBasedElement",
			  "parameters":{"filename":"/gws/nopw/j04/eustace/data/internal/climatology_covariates/DEM_global_0.25_0.25.nc",

Here users can define which physical quantities covary with the temperature field, along with the corresponding prior hyperparameters.

To include only the seasonal core into the climatology model, just substitute the json descriptor path with :code:`None`

Catalogue Generation
--------------------

To generate the catalogue of inputs and operations, you first need to create a directory where the analysis results will be stored.

You should look for storage spaces large enough to contain the output data. For example, you could create a directory called :code:`/work/scratch/${USERNAME}/advanced_standard`

Once an output directory has been created, execute the following command

.. code-block:: console

    $python2.7 -m eumopps.catalogue.build ${JSON_PATH}/advstd.json ${OUT_DIR}/catalogue.nc --pathdefault ${OUT_DIR}

where :code:`${JSON_PATH}` is the path to the folder containing :code:`advstd_inputs.json` and :code:`${OUT_DIR}` is the output folder you have previously created.

This command will build a catalogue containing information about all the inputs needed for the analysis system, for time periods specified by users.

If the catalogue creation has been successfull, you should find a netCDF file, called :code:`catalogue.nc` inside :code:`${OUT_DIR}`.

.. code-block:: console

    $ls ${OUT_DIR}/*nc
    catalogue.nc

If that is the case, then execute the following command

.. code-block:: console

    $python2.7 -m eumopps.catalogue.build --update ${JSON_PATH}/advstd_inputs.json ${OUT_DIR}/catalogue.nc --pathdefault ${OUT_DIR}

which will update the catalogue by adding information about the operations :ref:`EUMOPPS <eumopps>` will execute.

If you don't manage to produce a catalogue with the first command, then you should probably have a look at the changes you applied to the json descriptors.


Running the analyis
-------------------

Commands hierarchy
------------------

To run the analysis in the correct way, the following hierarchy of procedures should be respected:

1. **Climatology model**

   a. climatology input creation;
   b. climatology model solution.

2. **Large scale model**

  a. large scale input creation: this will use the climatology model solution to condition the large scale model;
  b. large scale input solution.

3. **Local input creation** and **model solution**: this will use the solutions from the higher level models to condition
   the local scale model and then solve it.

4. **Output production**: this will produce the final output onto a rectilinear grid, by combining the solutions obtained from the three temperature models mentioned above.

Commands execution
------------------

All these operations are performed by using the bash scripts contained in :code:`${HOME}/trunk/system/src/roses/fullstace/`

.. code-block:: console

   $ls ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/*sh
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/lsf_commandrun.sh
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/lsf_commandrun_solver.sh
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/runlsf.sh
   ${HOME}/trunk/system/src/roses/fullstace/advanced_standard/runlsf_solver.sh

More precisely, you will just need the last two scripts in the above list.

To run the analysis, move into :code:`${HOME}/trunk/system/src/roses/fullstace/advanced_standard` and execute the commands described below.

1. **Climatology model**

   a. climatology model: input creation
   
   .. code-block:: console

      $./runlsf.sh ${OUT_DIR}/catalogue.nc climatology_input
      Running: "./lsf_commandrun.sh" for EUMOPPS module "climatology_input" with 1 batches of 366                                                              
      Job <XXXXXX> is submitted to queue <short-serial>

   This will create a new folder into :code:`${OUT_DIR}`, called :code:`measurement_climatology`. The latter will contain a collection of directories labelled by the years covered by the analysis.

   Inside each folder you will find a set of :code:`pickle` files, one for all the day processed within a given year. These files store information about measurement precision and measurement vector.

   After having launched the production of climatology inputs, you can submit the other operations in queue, you will just need to copy the job ID :code:`XXXXXX`
   
   b. climatology model: solution

   .. code-block:: console

      $ ./runlsf_solve.sh ${OUT_DIR}/catalogue.nc climatology_solve "done(XXXXXX)"
      Running: "./lsf_commandrun_solve.sh" for EUMOPPS module "climatology_solve" with 1 batches of 1 
      Dependencies: "done(XXXXXX)"
      Job <YYYYYY> is submitted to queue <par-single>

   The script will generate a folder called :code:`climatology_solution`, with a :code:`pickle` file storing the posterior precision and the MAP estimates
   of the state parameters used for the climatology model.

   The third argument assigned to :code:`runlsf_solve.sh` will force the :code:`lsf` scheduler to wait the conclusion of job :code:`XXXXXX` before submitting the 
   computation of the climatology model solution.

   In this way, you can submit the execution of all the operations at once, being sure that the hierarchy of commands is respected.
   
2. **Large scale model**

  a. large scale model: input creation

  .. code-block:: console

      $./runlsf.sh ${OUT_DIR}/catalogue.nc large_scale_input "done(YYYYYY)"
      Running: "./lsf_commandrun.sh" for EUMOPPS module "large_scale_input" with 1 batches of 366                                                              
      Dependencies: "done(YYYYYY)"
      Job <WWWWWW> is submitted to queue <short-serial>

  b. large scale model: solution

  .. code-block:: console

      $ ./runlsf_solve.sh ${OUT_DIR}/catalogue.nc large_scale_solve "done(WWWWWW)"
      Running: "./lsf_commandrun_solve.sh" for EUMOPPS module "large_scale_solve" with 1 batches of 1 
      Dependencies: "done(WWWWW)"
      Job <ZZZZZZ> is submitted to queue <par-single>

  as for the climatology model case, two new directories, with :code:`pickle` files, will be created, called :code:`measurement_large_scale`, :code:`large_scale_solution`.

3. **Local input creation** and **model solution**

  .. code-block:: console

      $./runlsf.sh ${OUT_DIR}/catalogue.nc local_input_and_solve "done(ZZZZZZ)"
      Running: "./lsf_commandrun.sh" for EUMOPPS module "local_input_and_solve" with 1 batches of 366                                                              
      Dependencies: "done(ZZZZZZ)"
      Job <AAAAAA> is submitted to queue <short-serial>

  at the end of the computation, results will be stored into a :code:`solution_local` folder.

4. **Output production**

  .. code-block:: console

      $./runlsf.sh ${OUT_DIR}/catalogue.nc output_grid "done(AAAAAA)"
      Running: "./lsf_commandrun.sh" for EUMOPPS module "output_grid" with 1 batches of 366                                                              
      Dependencies: "done(AAAAAA)"
      Job <BBBBBB> is submitted to queue <short-serial>

  at the end of the computation, results will be stored into a :code:`eustace_example_infilled` folder.

  Here you will find a collection of directories, one for each year, with the results of the infilled analysis.

  .. code-block:: console

      $ls ${OUT_DIR}/eustace_example_infilled
      Y1 Y2 ... YN
      $ls ${OUT_DIR}/eustace_example_infilled/2006
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061226.nc
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061227.nc
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061228.nc
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061229.nc
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061230.nc
      ${OUT_DIR}/eustace_example_infilled/2006/eustace_example_infilled_20061231.nc

      

