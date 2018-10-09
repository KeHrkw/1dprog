MODULE CONSTANTS
  implicit none

  integer :: igauge ! 1=Velocity,2=Length,3=coeff
  integer :: Nx
  integer :: Nk
  integer :: Nb, Niter, Nscf

  real(8) :: E0_in, lambda_in, tau_in         !input parameter
  integer :: Nbt, Ne, Nexp
  real(8) ::dt, de

  real(8),parameter :: pi=acos(-1d0)
  real(8),parameter :: width=8.d0, a0=1.d0, V0=-0.37d0   !potential parameter
  complex(8),parameter :: zI=(0d0,1d0)
  !real(8),parameter :: Tad=500d0
  real(8),parameter :: t_au=0.0242d0 !fs
  real(8),parameter :: E0_au=5.14d0*1.d+2 !V/nm
  real(8),parameter :: a0_au=5.292*1.d-2 !nm
  real(8),parameter :: alpha_au=1d0/137d0
  real(8),parameter :: c_light=137.03953250d0
  real(8),parameter :: En_au=27.21d0 !eV


  real(8) :: dx, dk

  real(8) :: E0,omega,tau

  real(8),dimension(:),allocatable :: pot
  complex(8),dimension(:,:,:),allocatable :: phase

  integer,parameter :: Nd=4
  real(8),dimension(-Nd:Nd) :: lap,nab
  integer,dimension(:),allocatable :: modx
  integer :: NUMBER_THREADS
END MODULE CONSTANTS

MODULE WAVE_FUNC
  use CONSTANTS
  implicit none
  complex(8),dimension(:,:),allocatable :: zu_in
  complex(8),dimension(:,:),allocatable :: hzu, czu
  real(8),dimension(:),allocatable :: k
  complex(8),dimension(:,:,:),allocatable :: u
  complex(8),dimension(:,:,:) ,allocatable:: zu
  real(8) :: k_in, phase_in
  real(8),dimension(:),allocatable :: hav_base
END MODULE WAVE_FUNC

MODULE TD_CALC
  use CONSTANTS
  implicit none
  integer :: Nt
  real(8),dimension(:),allocatable  :: mask
  real(8),dimension(:,:),allocatable  :: cur, hav
  real(8),dimension(:,:),allocatable  ::  norm
  real(8),dimension(:),allocatable  :: Et, At

END MODULE TD_CALC
MODULE FD_K
  use CONSTANTS
  implicit none
  integer,parameter :: Nd_k=5
  real(8),dimension(-Nd_k:Nd_k) :: nab_k
  complex(8),dimension(:,:,:),allocatable :: dk_zu
  complex(8),dimension(:,:,:),allocatable :: zu_in_L
END MODULE

MODULE COEFF
  use CONSTANTS
  use FD_K
  implicit none
  integer :: Norb
  complex(8),dimension(:,:,:),allocatable :: Cnm
  complex(8),dimension(:,:),allocatable :: Cnm_in
  complex(8),dimension(:,:),allocatable :: HCnm
  complex(8),dimension(:,:,:),allocatable :: udku
  complex(8),dimension(:,:,:),allocatable :: udxu
  real(8),dimension(:,:),allocatable :: eps
END MODULE
