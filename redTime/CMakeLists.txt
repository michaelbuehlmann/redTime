set(REDTIME_SOURCES redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h)

add_executable(redTime ${REDTIME_SOURCES})
target_link_libraries(redTime gsl gslcblas)