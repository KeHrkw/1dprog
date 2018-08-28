SUBROUTINE Current_Velocity_operation()
  !Current operator
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  integer :: ix
  !complex(kind(0d0)) :: dif
# define IDX(dt) modx(ix+(dt)+Nx)
  complex(8) :: w
  !$omp threadprivate(czu_V,zu_in_V,k_in)

  !zu_in_V(-1)  =zu_in_V(Nx-1)
  !zu_in_V(Nx)  =zu_in_V(0)
  !zu_in_V(-2)  =zu_in_V(Nx-2)
  !zu_in_V(Nx+1)=zu_in_V(1)
  !zu_in_V(-3)  =zu_in_V(Nx-3)
  !zu_in_V(Nx+2)=zu_in_V(2)
  !zu_in_V(-4)  =zu_in_V(Nx-4)
  !zu_in_V(Nx+3)=zu_in_V(3)

  do ix = 0 ,Nx-1
    !dif = (zu_in_V(ix+1)-zu_in_V(ix-1))/(2.d0*dx)
    !dif = (-zu_in_V(ix+2) +8.d0*zu_in_V(ix+1) -8.d0*zu_in_V(ix-1)+zu_in_V(ix-2))/(12.d0*dx)
    !dif = (-224.d0*zu_in_V(ix+1)+224.d0*zu_in_V(ix-1)+56.d0*zu_in_V(ix+2)-56.d0*zu_in_V(ix-2) &
    !      & -32.d0/3.d0*zu_in_V(ix+3)+32.d0/3.d0*zu_in_V(ix-3)+zu_in_V(ix+4)-zu_in_V(ix-4))/(-280.d0*dx)

    w=(nab(1)*(zu_in_V(IDX(1))-zu_in_V(IDX(-1))) &
    & +nab(2)*(zu_in_V(IDX(2))-zu_in_V(IDX(-2))) &
    & +nab(3)*(zu_in_V(IDX(3))-zu_in_V(IDX(-3))) &
    & +nab(4)*(zu_in_V(IDX(4))-zu_in_V(IDX(-4))))

    czu_V(ix)= (-zI)*w + k_in*zu_in_V(ix) + At_in*zu_in_V(ix)
  end do

END SUBROUTINE
