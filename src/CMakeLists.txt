set(REDTIME_SOURCES redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h)

add_executable(redTime ${REDTIME_SOURCES})
target_link_libraries(redTime PRIVATE GSL::gsl GSL::gslcblas)
if(OpenMP_CXX_FOUND)
    target_link_libraries(redTime PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(redTime PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/")