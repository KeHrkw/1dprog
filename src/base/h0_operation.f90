SUBROUTINE h0_operation()
  !Hamiltonian of Grand-state
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  integer :: ix
  !complex(8) :: dif, dif2
  !complex(8),intent(in)  :: E(0:Nx-1)
  !complex(8),intent(out) :: F(0:Nx-1)
# define IDX(dt) modx(ix+(dt)+Nx)
  complex(8) :: v,w

  !u_in(-1)=u_in(Nx-1)
  !u_in(Nx)=u_in(0)
  !u_in(-2)=u_in(Nx-2)
  !u_in(Nx+1)=u_in(1)
  !u_in(-3)=u_in(Nx-3)
  !u_in(Nx+2)=u_in(2)
  !u_in(-4)=u_in(Nx-4)z
  !u_in(Nx+3)=u_in(3)
  do ix = 0, Nx-1, 1
    !dif = (u_in(ix+1)-u_in(ix-1))/(2.d0*dx)
    !dif = (-u_in(ix+2) +8.d0*u_in(ix+1) -8.d0*u_in(ix-1)+u_in(ix-2))/(12.d0*dx)
    !dif = (u_in(ix+3)-9.d0*u_in(ix+2)+45.d0*u_in(ix+1)-45.d0*u_in(ix-1)+9.d0*u_in(ix-2)-u_in(ix-3))/(60.d0*dx)
    !dif = (-224.d0*u_in(IDX(1))+224.d0*u_in(IDX(-1)) &
    !      & +56.d0*u_in(IDX(2))-56.d0*u_in(IDX(-2)) &
    !      & -32.d0/3.d0*u_in(IDX(3))+32.d0/3.d0*u_in(IDX(-3)) &
    !      &       +u_in(IDX(4))       -u_in(IDX(-4)))/(-280.d0*dx)

    !dif2 = (u_in(ix+1)-2.d0*u_in(ix)+u_in(ix-1))/(dx*dx)
    !dif2 =(-u_in(ix+2)+16.d0*u_in(ix+1)+16.d0*u_in(ix-1)-u_in(ix-2)-30.d0*u_in(ix))/(12.d0*dx*dx)
    !dif2 = (u_in(ix+3)+u_in(ix-3)-13.5d0*u_in(ix+2)-13.5d0*u_in(ix-2)+135.d0*u_in(ix+1)+135.d0*u_in(ix-1)-245.d0*u_in(ix))/(90.d0*dx*dx)
    !dif2 = (-896.d0*u_in(IDX(1))-896.d0*u_in(IDX(-1)) &
    !      & +112.d0*u_in(IDX(2))+112.d0*u_in(IDX(-2)) &
    !      & -128.d0/9.d0*u_in(IDX(3))-128.d0/9.d0*u_in(IDX(-3)) &
    !      &        +u_in(IDX(4))+u_in(IDX(-4)) &
    !      &    +2.d0*7175.d0/9.d0*u_in(ix))/(-560.d0*dx*dx)

    v=(lap(1)*(u_in(IDX(1))+u_in(IDX(-1))) &
    & +lap(2)*(u_in(IDX(2))+u_in(IDX(-2))) &
    & +lap(3)*(u_in(IDX(3))+u_in(IDX(-3))) &
    & +lap(4)*(u_in(IDX(4))+u_in(IDX(-4))))
    w=(nab(1)*(u_in(IDX(1))-u_in(IDX(-1))) &
    & +nab(2)*(u_in(IDX(2))-u_in(IDX(-2))) &
    & +nab(3)*(u_in(IDX(3))-u_in(IDX(-3))) &
    & +nab(4)*(u_in(IDX(4))-u_in(IDX(-4))))

    !hu(ix)=0.5d0*(-dif2 - 2.d0*zI*(k_in)*dif + (k_in)*(k_in)*u_in(ix))+pot(ix)*u_in(ix)
    hu(ix) = pot(ix)*u_in(ix) + 0.5d0*((k_in)*(k_in)-lap(0))*u_in(ix) - 0.5d0 * v - zI * k_in * w
  end do


END SUBROUTINE
