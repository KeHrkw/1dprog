SUBROUTINE zh_Velocity_operation(kAc,A,B)
  use CONSTANTS, only: Nx,nab,lap,zI,modx,pot
  implicit none
  integer :: ix
  !complex(kind(0d0)) :: dif, dif2
# define IDX(dt) modx(ix+(dt)+Nx)
  complex(8),dimension(0:Nx-1),intent(in)  :: A
  complex(8),dimension(0:Nx-1),intent(out) :: B
  real(8),intent(in)    :: kAc
  complex(8) :: v,w



  do ix = 0, Nx-1, 1

    v=(lap(1)*(A(IDX(1))+A(IDX(-1))) &
    & +lap(2)*(A(IDX(2))+A(IDX(-2))) &
    & +lap(3)*(A(IDX(3))+A(IDX(-3))) &
    & +lap(4)*(A(IDX(4))+A(IDX(-4))))
    w=(nab(1)*(A(IDX(1))-A(IDX(-1))) &
    & +nab(2)*(A(IDX(2))-A(IDX(-2))) &
    & +nab(3)*(A(IDX(3))-A(IDX(-3))) &
    & +nab(4)*(A(IDX(4))-A(IDX(-4))))

    B(ix) = pot(ix)*A(ix) + 0.5d0*((kAc)*(kAc)-lap(0))*A(ix) - 0.5d0 * v - zI * (kAc) * w
  end do

END SUBROUTINE
