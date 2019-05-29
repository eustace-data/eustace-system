Working with the C++ API
========================

The C++ API is built using the `CMake`_ tool which is part of the installed prequisites.
The simplest way to build new programs which link to the API is to use CMake for those too,
and cross-reference the EUSTACE library.

A first example
---------------

Create a directory ``eustace_tryout`` to contain example code:

.. code-block:: bash

    mkdir eustace_tryout
    cd eustace_tryout

And in this directory create a ``CMakeLists.txt`` and ``eustace_tryout.cpp`` as follows:

**CMakeLists.txt**

.. code-block:: bash

    cmake_minimum_required (VERSION 2.8)

    # Use Python to get EUSTACE code location
    execute_process(COMMAND python -c "import eustaceconfig;print eustaceconfig.CODE_PATH" OUTPUT_VARIABLE EUSTACE_CODE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Cross-reference build location
    include_directories("${EUSTACE_CODE_PATH}/cpp")
    link_directories("${EUSTACE_CODE_PATH}/build/eustace")

    # Define the executable
    add_executable(eustace_tryout eustace_tryout.cpp)

    # Link to the EUSTACE API
    target_link_libraries(eustace_tryout eustace)


**eustace_tryout.cpp**

.. code-block:: C++

    #include <eustace/timeutils/timebase.h>
    #include <eustace/definitions.h>
    #include <iostream>

    int main(int argc, char* argv[])
    {
      std::string datestring("20170222");
      EUSTACE::CalendarDay day(datestring);
      std::cout << 
	"Day number for date " << datestring << 
	" is " << EUSTACE::TimeBaseDays(EUSTACE::EPOCH).Number(day) << std::endl;
    }

Now make a ``build`` subdirectory and build it:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make

If all has worked as it should, it should now be possible to run the program ``./eustace_tryout``:

.. code-block:: text

    $ ./eustace_tryout 
    Day number for date 20170222 is 61048

.. _CMake: https://cmake.org/
