set(REDTIME_SOURCES redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h)

add_executable(redTime ${REDTIME_SOURCES})
target_link_libraries(redTime PRIVATE GSL::gsl GSL::gslcblas)
if(OpenMP_CXX_FOUND)
    target_link_libraries(redTime PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(redTime PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/")

add_executable(convertPk convert_pk.c)
target_link_libraries(convertPk PRIVATE GSL::gsl GSL::gslcblas m)
set_target_properties(convertPk PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/")
