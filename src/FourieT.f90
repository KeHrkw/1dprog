SUBROUTINE FourieT()
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  complex(kind(0d0)),dimension(:) :: ecur(Ne)
  real(8) :: ee, tt, ss
  integer :: ie, it

  write(*,*)
  write(*,*) "Fourie Transform"
  write(*,*)
  write(*,*)"de : ", de
  write(*,*)"Ne : ", Ne

  ecur(:)=0.d0
  !$omp parallel do private(ie,ee,it,tt,ss)
  do ie = 1, Ne, 1
    ee = ie*de
      do it = 1, Nt, 1
        tt = it*dt
        ss = sum(cur(1:Nbt,it))/real(Nbt)
        ecur(ie) = ecur(ie) + exp(zI*tt*ee)* ss *smoothing_t(tt)
      end do
      ecur(ie) = ecur(ie)*dt/(zI*ee)
  end do
  !$omp end parallel do

  open(9,file='./fts.data')
    do ie = 1, Ne, 1
      ee = ie*de
      write(9,'(<1>i,<3>e)') ie,ee,ee/omega,abs(ecur(ie))**2
    end do
  close(9)
contains
  Function smoothing_t(tt)
    implicit none
    real(8),intent(IN) :: tt
    real(8) :: smoothing_t

    smoothing_t=1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3
    !if(abs(tt)<0.5d0*pulse_tw1) then
    !  smoothing_t=cos(pi*tt/pulse_tw1)**2
    !end if

    return
  End Function smoothing_t
END SUBROUTINE
