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
    & E0_in, &
    & lambda_in, &
    & tau_in, &
    & Nbt, &
    & Nexp, &
    & dt, &
    & de, &
    & Ne

!! == default for &paramet
  igauge=1
  NX=15
  Nk=21
  Nb=3
  Niter=120
  Nscf=5
  E0_in=1.65d0
  lambda_in=3200d0
  tau_in=42.66d0
  Nbt=2
  Nexp=4
  dt=0.011d0
  de=0.02d0
  Ne=500

  open(10,file='input.data')
  read(10,nml=group)
  write(*,nml=group)
  close(10)

  dx=width/real(Nx)
  dk=pi/width/(real(Nk))*2.d0
  E0=E0_in/E0_au
  omega=a0_au/lambda_in*2d0*pi*c_light
  tau=tau_in/t_au

  Nt=int(pi*tau/dt)

  Norb=Nb
end subroutine
