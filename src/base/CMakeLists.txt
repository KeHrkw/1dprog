set(myapp_sources
  base_plot.f90
  CG.f90
  h0_operation.f90
  Imaginary_time.f90
  Phase_matching.f90
)

add_library(BASE_LIB ${myapp_sources} )
