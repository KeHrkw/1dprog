SUBROUTINE FourieT()
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  complex(kind(0d0)),dimension(:,:) :: ecur(Nbt,Ne)
  real(8) :: peak
  real(8) :: ee, tt, ss
  integer :: ie, dw, ib, jb, it

  write(*,*)
  write(*,*) "Fourie Transform"
  write(*,*)
  write(*,*)"de : ", de
  write(*,*)"Ne : ", Ne

  ecur(:,:)=0.d0
  !$omp parallel do private(ie,ee,ib,it,tt,ss)
  do ie = 1, Ne, 1
    ee = ie*de
    do ib = 1, Nbt, 1
      do it = 1, Nt, 1
        tt = it*dt
        ss = sum(cur(1:ib,it))
        ecur(ib,ie) = ecur(ib,ie) + exp(zI*tt*ee)* ss *mask(it)
      end do
      ecur(ib,ie) = ecur(ib,ie)*dt/(zI*ee)
    end do
  end do
  !$omp end parallel do

  open(9,file='./output/fts.dat')
  do ib = 1, Nbt, 1
    do ie = 1, Ne, 1
      ee = ie*de
      write(9,'(<2>i,<3>e)') ib,ie,ee,ee*En_au,abs(ecur(ib,ie))**2/ib
    end do
    write(9,*)""
    write(9,*)""
  end do
  close(9)

  dw=int(omega/de)
  write(*,*)"dw : ", dw
  write(*,*)"Ne/dw : ", Ne/dw
  open(8,file='./output/integ.dat')
  do ib = 1, Nbt, 1
    do jb = 1, Ne/dw-1, 2
      peak=0.d0
      do ie = (jb-1)*dw+1, (jb+1)*dw, 1
        peak = peak +abs(ecur(ib,ie))**2/ib*de
      end do
      write(8,'(<2>i,<2>e)') ib,jb,jb*dw*de,peak
    end do
  write(8,*) ""
  write(8,*) ""

  end do
  close(8)

END SUBROUTINE
