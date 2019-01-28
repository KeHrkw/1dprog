subroutine read_input_file
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  use COEFF
  use FD_K
  implicit none
  namelist / group / &
    & igauge, &
    & Nx, &
    & Nk, &
    & Nb, &
    & Niter, &
    & Nscf, &
    & rlaser_int_wcm2, &
    & lambda_in, &
    & tau_in, &
    & Nt, &
    & Nbt, &
    & Nexp, &
    & dt, &
    & de, &
    & Ne
  namelist / potential / &
    & ipot, &
    & width, &
    & a0, &
    & V0

!! == default for &paramet
  igauge=1
  NX=15
  Nk=21
  Nb=3
  Niter=120
  Nscf=5
  rlaser_int_wcm2=0d0
  lambda_in=3200d0
  tau_in=42.66d0
  Nt=3000
  Nbt=2
  Nexp=4
  dt=0.011d0
  de=0.02d0
  Ne=500


  ipot=1
  width=8d0
  a0=1d0
  V0=-7.d0

  open(10,file='input.data')
  read(10,nml=group)
  write(*,nml=group)
  read(10,nml=potential)
  write(*,nml=potential)
  close(10)

  dx=width/real(Nx)
  dk=pi/width/(real(Nk))*2.d0
  !E0=E0_in/E0_au
  E0=5.338d-9*sqrt(rlaser_int_wcm2)
  omega=a0_au/lambda_in*2d0*pi*c_light
  tau=tau_in
  !tau=tau_in/t_au

  !Nt=int(pi*tau/dt)

  Norb=Nb
end subroutine
