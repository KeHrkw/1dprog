set(myapp_sources
  # Current_Length_operation.f90
  Gauge_Length.f90
  # zh_Length_operation.f90
)
add_library(LENGTH_LIB ${myapp_sources} )
message("LENGTH_LIB = ${LENGTH_LIB}")
