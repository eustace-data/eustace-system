# Windows configuration parameters
if (MSVC)
  option(EUSTACE_MSVC_RUNTIME_STATIC "Use static (non-DLL) version of C runtime libraries" ON)
endif(MSVC)

# Setup paths on cems
option(ON_CEMS "Running on CEMS" ON)

if (ON_CEMS)
  include_directories("/apps/libs/netCDF/intel17/4.6.1/include")
  link_directories("/apps/libs/netCDF/intel17/4.6.1/lib/")
endif()

# Optionally use static libraries on windows
if (EUSTACE_MSVC_RUNTIME_STATIC)
  foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO 
	CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
    if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif(${flag_var} MATCHES "/MD")
  endforeach(flag_var)
endif(EUSTACE_MSVC_RUNTIME_STATIC)

add_library(eustace
  definitions.cpp
  timeutils/timebase.cpp
  analysis/observationsource.cpp
  analysis/inputmanager.cpp
  analysis/fileio/rawbinary.cpp
  analysis/fileio/locationpointers.cpp
  analysis/fileio/mesh.cpp
  analysis/break_points/break_points.cpp
  analysis/break_points/insitu_land_break_points.cpp
  )
add_subdirectory(test)
add_subdirectory(examples)
