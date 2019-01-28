set(myapp_sources
  Gauge_Length.f90
)
add_library(LENGTH_LIB ${myapp_sources} )
message("LENGTH_LIB = ${LENGTH_LIB}")
