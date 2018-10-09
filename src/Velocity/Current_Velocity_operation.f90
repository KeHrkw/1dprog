SUBROUTINE Current_Velocity_operation(kAc,A,B)
  !Current operator
  use CONSTANTS, only: Nx,nab,zI,modx
  implicit none
  integer :: ix
  !complex(kind(0d0)) :: dif
# define IDX(dt) modx(ix+(dt)+Nx)
  complex(8),dimension(0:Nx-1),intent(in)  :: A
  complex(8),dimension(0:Nx-1),intent(out) :: B
  real(8),intent(in) :: kAc
  complex(8) :: w

  do ix = 0 ,Nx-1

    w=(nab(1)*(A(IDX(1))-A(IDX(-1))) &
    & +nab(2)*(A(IDX(2))-A(IDX(-2))) &
    & +nab(3)*(A(IDX(3))-A(IDX(-3))) &
    & +nab(4)*(A(IDX(4))-A(IDX(-4))))

    B(ix)= (-zI)*w + kAc*A(ix)
  end do

END SUBROUTINE
