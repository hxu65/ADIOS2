#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

add_test(NAME HeatTransfer.SST.BP.Mx1
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}
    ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:adios2_simulations_heatTransferWrite>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_sst_bp.xml
        WriteSSTBPMx1.bp 2 2 10 10 10 10 SST
    :
    ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:adios2_simulations_heatTransferRead>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_sst_bp.xml
        WriteSSTBPMx1.bp ReadSSTBPMx1.bp 1 1 SST
)
set_tests_properties(HeatTransfer.SST.BP.Mx1 PROPERTIES PROCESSORS 5)

add_test(NAME HeatTransfer.SST.BP.Mx1.Dump
  WORKING_DIRECTORY SSTBPMx1
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-d 
    -DINPUT_FILE=ReadSSTBPMx1.bp
    -DOUTPUT_FILE=DumpSSTBPMx1.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME HeatTransfer.SST.BP.Mx1.Validate
  WORKING_DIRECTORY SSTBPMx1
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/HeatTransfer.Dump.txt
    DumpSSTBPMx1.txt
)

SetupTestPipeline(HeatTransfer.SST.BP.Mx1 ";Dump;Validate" TRUE)
