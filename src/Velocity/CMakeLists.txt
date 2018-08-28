set(myapp_sources
  Current_Velocity_operation.f90
  Gauge_Velocity.f90
  zh_Velocity_operation.f90
)
add_library(VELOCITY_LIB ${myapp_sources} )
