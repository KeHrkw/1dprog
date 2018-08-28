MODULE CONSTANTS
  implicit none

  !integer,parameter :: LNk_part=5
  !integer,parameter :: LNk_1=100, LNk_2=5
  !real(8) :: dk_1, dk_2

  integer,parameter :: igauge=1! 1=Velocity,2=Length,3=coeff
  integer,parameter :: Nx=30, LNk=20, RNk=20
  integer,parameter :: Nk=LNk+RNk+1
  integer,parameter :: Nb=3, Niter=120, Nscf=5

  real(8),parameter :: E0_in=1.65d0, lambda_in=3200d0, tau_in=99.66d0          !input parameter
  integer,parameter :: Nbt=2, Ne=5000, Nexp=4
  real(8),parameter :: pi=acos(-1d0)
  real(8),parameter :: width=8.d0, a0=1.d0, V0=-0.37d0, dt=0.011d0, de=0.002d0
  complex(kind(0d0)),parameter :: zI=(0d0,1d0)
  !real(8),parameter :: Tad=500d0
  real(8),parameter :: dx=width/real(Nx), dk=pi/width/(real(Nk))*2.d0

  real(8),parameter :: t_au=0.0242d0 !fs
  real(8),parameter :: E0_au=5.14d0*1.d+2 !V/nm
  real(8),parameter :: a0_au=5.292*1.d-2 !nm
  real(8),parameter :: alpha_au=1d0/137d0
  real(8),parameter :: En_au=27.21d0 !eV

  real(8),parameter :: E0=E0_in/E0_au
  real(8),parameter :: omega=a0_au/alpha_au/lambda_in
  real(8),parameter :: tau=tau_in/t_au

  complex(kind(0d0)),dimension(0:Nx-1,1:Nb,0:1) :: phase

  integer,parameter :: Nd=4
  real(8),dimension(-Nd:Nd) :: lap,nab
  integer,dimension(:),allocatable :: modx
  !$ integer :: nthreads
END MODULE CONSTANTS

MODULE WAVE_FUNC
  use CONSTANTS
  implicit none
  real(8),dimension(:) :: pot(0:Nx-1)
  real(8),dimension(:) :: k(-LNk:RNk)
  complex(kind(0d0)),dimension(:,:,:) :: u(0:Nx-1,-LNk:RNk,Nb)
  complex(kind(0d0)),dimension(:,:,:) :: zu(0:Nx-1,-LNk:RNk,Nbt)
  complex(kind(0d0)),dimension(:) :: u_in(0:Nx-1)
  complex(kind(0d0)),dimension(:) ::   hu(0:Nx-1)
  real(8) :: k_in, phase_in
  complex(kind(0d0)),dimension(:) :: zu_in_V(0:Nx-1)
  complex(kind(0d0)),dimension(:) :: hzu_V(0:Nx-1), czu_V(0:Nx-1)
  !complex(kind(0d0)),dimension(:,:) :: zu_in_L(-4:Nx+3,-LNk-5:RNk+5)
  !complex(kind(0d0)),dimension(:,:) :: hzu_L(0:Nx-1,-LNk:RNk), czu_L(0:Nx-1,-LNk:RNk)
  real(8),dimension(:) :: hav_base(Nb)
END MODULE WAVE_FUNC

MODULE TD_CALC
  use CONSTANTS
  implicit none
  integer,parameter :: Nt=int(pi*tau/dt)
  real(8),dimension(0:Nt) :: mask
  real(8),dimension(1:Nbt,0:Nt) :: cur, hav
  real(8),dimension(1:Nbt) ::  norm
  real(8),dimension(0:Nt) :: Et, At
  real(8) :: At_in, Et_in
  integer :: it

END MODULE TD_CALC
MODULE FD_K
  use CONSTANTS
  implicit none
  integer,parameter :: Nd_k=5
  real(8),dimension(-Nd_k:Nd_k) :: nab_k
  complex(8),dimension(:,:) :: dk_zu(0:Nx-1,-LNk:RNk)
  complex(kind(0d0)),dimension(:,:) :: zu_in_L(0:Nx-1,-LNk-Nd_k:RNk+Nd_k)
END MODULE

MODULE COEFF
  use CONSTANTS
  use FD_K
  implicit none
  integer,parameter :: Norb=Nb
  complex(kind(0d0)),dimension(:,:,:) :: Cnm(-LNk:RNk,Nbt,Norb)
  complex(kind(0d0)),dimension(:,:,:) :: Cnm_in(-LNk-Nd_k:RNk+Nd_k,Norb)
  complex(kind(0d0)),dimension(:,:,:) :: HCnm(-LNk:RNk,Norb)
  complex(kind(0d0)),dimension(:,:,:) :: udku(-LNk:RNk,Norb,Norb)
  complex(kind(0d0)),dimension(:,:,:) :: udxu(-LNk:RNk,Norb,Norb)
  real(8),dimension(:,:) :: eps(-LNk:RNk,Nb)
END MODULE
