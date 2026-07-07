#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "adios2::perfstubs" for configuration "Release"
set_property(TARGET adios2::perfstubs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::perfstubs PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libadios2_perfstubs.so.2.11.0"
  IMPORTED_SONAME_RELEASE "libadios2_perfstubs.so.2.11"
  )

list(APPEND _IMPORT_CHECK_TARGETS adios2::perfstubs )
list(APPEND _IMPORT_CHECK_FILES_FOR_adios2::perfstubs "${_IMPORT_PREFIX}/lib/libadios2_perfstubs.so.2.11.0" )

# Import target "adios2::core" for configuration "Release"
set_property(TARGET adios2::core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::core PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "adios2::core_derived;Blosc2::blosc2_shared;zfp::zfp;mgard::mgard"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libadios2_core.so.2.11.0"
  IMPORTED_SONAME_RELEASE "libadios2_core.so.2.11"
  )

list(APPEND _IMPORT_CHECK_TARGETS adios2::core )
list(APPEND _IMPORT_CHECK_FILES_FOR_adios2::core "${_IMPORT_PREFIX}/lib/libadios2_core.so.2.11.0" )

# Import target "adios2::core_mpi" for configuration "Release"
set_property(TARGET adios2::core_mpi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::core_mpi PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libadios2_core_mpi.so.2.11.0"
  IMPORTED_SONAME_RELEASE "libadios2_core_mpi.so.2.11"
  )

list(APPEND _IMPORT_CHECK_TARGETS adios2::core_mpi )
list(APPEND _IMPORT_CHECK_FILES_FOR_adios2::core_mpi "${_IMPORT_PREFIX}/lib/libadios2_core_mpi.so.2.11.0" )

# Import target "adios2::core_derived" for configuration "Release"
set_property(TARGET adios2::core_derived APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::core_derived PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libadios2_core_derived.so"
  IMPORTED_SONAME_RELEASE "libadios2_core_derived.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS adios2::core_derived )
list(APPEND _IMPORT_CHECK_FILES_FOR_adios2::core_derived "${_IMPORT_PREFIX}/lib/libadios2_core_derived.so" )

# Import target "adios2::ParaViewADIOSInSituEngine" for configuration "Release"
set_property(TARGET adios2::ParaViewADIOSInSituEngine APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::ParaViewADIOSInSituEngine PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libParaViewADIOSInSituEngine.so"
  IMPORTED_NO_SONAME_RELEASE "TRUE"
  )

list(APPEND _IMPORT_CHECK_TARGETS adios2::ParaViewADIOSInSituEngine )
list(APPEND _IMPORT_CHECK_FILES_FOR_adios2::ParaViewADIOSInSituEngine "${_IMPORT_PREFIX}/lib/libParaViewADIOSInSituEngine.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
