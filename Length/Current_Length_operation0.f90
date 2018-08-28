SUBROUTINE Current_Length_operation()
  !Current operator
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  !$omp threadprivate(zu_in_L,czu_L)

  integer :: ix, ik
  complex(kind(0d0)) :: dif

  zu_in_L(-1,:)  =zu_in_L(Nx-1,:)
  zu_in_L(Nx,:)  =zu_in_L(0,:)
  zu_in_L(-2,:)  =zu_in_L(Nx-2,:)
  zu_in_L(Nx+1,:)=zu_in_L(1,:)
  zu_in_L(-3,:)  =zu_in_L(Nx-3,:)
  zu_in_L(Nx+2,:)=zu_in_L(2,:)
  zu_in_L(-4,:)  =zu_in_L(Nx-4,:)
  zu_in_L(Nx+3,:)=zu_in_L(3,:)
  do ik = -LNk, RNk, 1
    do ix = 0 ,Nx-1
      !dif=(zu_in_L(ix+1,ik)-zu_in_L(ix-1,ik))/(2.d0*dx)
      !dif =(-zu_in_L(ix+2,ik) +8.d0*zu_in_L(ix+1,ik) -8.d0*zu_in_L(ix-1,ik)+zu_in_L(ix-2,ik))/(12.d0*dx)
      dif = (-224.d0*zu_in_L(ix+1,ik)+224.d0*zu_in_L(ix-1,ik)+56.d0*zu_in_L(ix+2,ik)-56.d0*zu_in_L(ix-2,ik) &
            & -32.d0/3.d0*zu_in_L(ix+3,ik)+32.d0/3.d0*zu_in_L(ix-3,ik)+zu_in_L(ix+4,ik)-zu_in_L(ix-4,ik))/(-280.d0*dx)
      czu_L(ix,ik)= (-zI)*dif + k(ik)*zu_in_L(ix,ik)
    end do
  end do

END SUBROUTINE
