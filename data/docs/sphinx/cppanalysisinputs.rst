Analysis input data available from C++
======================================

Example usage is demonstrated in the example program :ref:`cpp_example_analysis_input`.
Observation input data, correlation ranges, and location information are obtained using the
:cpp:class:`AnalysisInputManager`.

At startup, the calling program should call the :cpp:member:`AnalysisInputManager::Specify` 
method to indicate the observations of interest. To retrieve observations at a fixed time, 
use :cpp:func:`AnalysisInputManager::RetrieveDay`,
or to retrieve observations for all time at a fixed set of locations, 
use :cpp:func:`AnalysisInputManager::RetrieveLocations`.

The results of a daily retrieval are expressed in a :cpp:class:`AnalysisInputRetrieval` object,
and the results of a location-based retrieval are expressed as a collection 
of :cpp:class:`AnalysisInputRetrieval` objects, in the form of :cpp:class:`AnalysisInputRetrievalCollection`.

Observations contained within an :cpp:class:`AnalysisInputRetrieval` can be accessed using the iterator class
:cpp:class:`AnalysisInputIterator`, which also provides facility to determine the location of the observations.

Example Programs
----------------

.. toctree::

   cpp_example_analysis_input.rst

Use of Location ID
------------------

Internally, data sources have one or more lookup tables which convert from location identifiers (location ID) to geographic (latitude, longitude) coordinates.
These lookup tables are loaded automatically as described in the documentation for :cpp:func:`AnalysisInputManager::RetrieveDay` and 
:cpp:func:`AnalysisInputManager::RetrieveLocations`.

:cpp:func:`AnalysisInputIterator::Location` automatically performs the lookup and returns geographic coordinates. The location ID number for any observation can be accessed directly using :cpp:func:`AnalysisInputIterator::Identifier` if required.  This may be useful for precomputing other information that depends on location.  A location ID is a 64-bit integer.

For land station and satellite-derived data, the mapping between location ID and geographic location does not change with time.
For these sources, IDs will be in the range 0 to (:cpp:func:`AnalysisInputManager::TotalLocations` - 1) inclusive.  To find the coordinates corresponding
to a given location ID, it suffices to load any day of data using :cpp:func:`AnalysisInputManager::RetrieveDay`
and use the :cpp:func:`AnalysisInputRetrieval::ComputeLocationLookup` method.

Ship data is a mobile data source.  For mobile data sources, there is a different location lookup table for each day.  IDs for a given day of ship data will be in the range
0 to (:cpp:func:`AnalysisInputManager::TotalDailyMobileLocations` - 1) inclusive.  Again the conversion to coordinates can be done using :cpp:func:`AnalysisInputManager::RetrieveDay` followed by :cpp:func:`AnalysisInputRetrieval::ComputeLocationLookup`, but in this case it is important that the location lookup for the day of interest is used.  Location IDs do not have the same meaning from one day to another for mobile data sources.


Class Reference
---------------

.. cpp:class:: AnalysisInputManager

  .. cpp:function:: AnalysisInputManager::AnalysisInputManager (const char* pathname)

      Construct an analysis input manager to access raw binary files on the specified path.

      No memory is allocated at this stage.

  .. cpp:function:: AnalysisInputManager::Specify(AnalysisInput& input)

      Should be called during a startup phase to tell the input manager which data sources are of
      interest in case preparations are required. At present this does nothing but record the list
      of requests.

      No data is loaded or memory allocated.

  .. cpp:function:: void AnalysisInputManager::RetrieveDay(AnalysisInputRetrieval& retrieval, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber)

      Retrieves all data for the given source, observable, and day number.
    
      Memory for observations is allocated in the retrieval object and observations read from disk.
    
      For data sources with mobile locations (e.g. ships), memory is allocated for the location coordinates,
      and these are read from disk.
    
      For data sources which are fixed for all time, memory allocation and disk-reading of location coordinates
      happens only during the first call to this method for a given data source and observable.
      Subsequent calls obtain a cached copy.  This cache is only de-allocated when the input manager is destroyed.
    
      Reusing the same retrieval object instance for further calls to this method will cause any 
      observations or mobile location information it previously held to be destroyed.
      Otherwise this information will be destroyed whenever the retrieval object is destroyed.
    
      Local correlation range information is read from disk the first time a data source and observable 
      is accessed.  After that a cache is used until the input manager is destroyed.

  .. cpp:function:: void AnalysisInputManager::RetrieveLocations(AnalysisInputRetrievalCollection& collection, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, std::vector<uint64_t>& location_id_subset)

      Retrieves all data for the given source, observable, and array of locations.
      This loads only data sources whose locations are fixed (i.e. not ship data).

      Memory allocation scheme is similar to :cpp:func:`AnalysisInputManager::RetrieveDay`.  
      The collection object contains multiple arrays of observation data, one array per location.
      These are allocated and read from disk when this method is called.
      They are destroyed whenever the collection object is destroyed.
      As with :cpp:func:`AnalysisInputManager::RetrieveDay`, the locations coordinates and correlation ranges are read the first 
      time that the source and observable are accessed, and cached after that.

  .. cpp:function:: size_t AnalysisInputManager::TotalLocations(enum AnalysisInput::sourcetype source) const

      Total locations for the specified source.  Applies to fixed sources only (i.e. not ships).

      This information is read from location coordinate file header to avoid time reading the data itself.

  .. cpp:function:: size_t AnalysisInputManager::TotalDailyMobileLocations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const

      Total locations for the specified source and observable on the given day.  Applies to mobile sources only (e.g. ships).

      This information is read from location coordinate file header for the given source/observable/day, to avoid time reading the data itself.

  .. cpp:function:: size_t AnalysisInputManager::TotalObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable)

      Total observations for given source and observable.

     Uses multiple calls to :cpp:func:`AnalysisInputManager::TotalDailyObservations` so accesses multiple file headers but does not read any observation data.

  .. cpp:function:: size_t AnalysisInputManager::TotalDailyObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber)

      Total observations for given source and observable on the given day.

      This information is read from observation data file header to avoid time reading the data itself.

