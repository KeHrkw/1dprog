SUBROUTINE zh_Velocity_operation()
  !Hamiltonian of Grand-state
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  integer :: ix
  !complex(kind(0d0)) :: dif, dif2
# define IDX(dt) modx(ix+(dt)+Nx)
  complex(8) :: v,w
  !$omp threadprivate(hzu_V,zu_in_V,k_in)

  !zu_in_V(-1)  =zu_in_V(Nx-1)
  !zu_in_V(Nx)  =zu_in_V(0)
  !zu_in_V(-2)  =zu_in_V(Nx-2)
  !zu_in_V(Nx+1)=zu_in_V(1)
  !zu_in_V(-3)  =zu_in_V(Nx-3)
  !zu_in_V(Nx+2)=zu_in_V(2)
  !zu_in_V(-4)  =zu_in_V(Nx-4)
  !zu_in_V(Nx+3)=zu_in_V(3)


  do ix = 0, Nx-1, 1
    !dif=(zu_in_V(ix+1)-zu_in_V(ix-1))/(2.d0*dx)
    !dif =(-zu_in_V(ix+2) +8.d0*zu_in_V(ix+1) -8.d0*zu_in_V(ix-1)+zu_in_V(ix-2))/(12.d0*dx)
    !dif = (-224.d0*zu_in_V(ix+1)+224.d0*zu_in_V(ix-1)+56.d0*zu_in_V(ix+2)-56.d0*zu_in_V(ix-2) &
    !      & -32.d0/3.d0*zu_in_V(ix+3)+32.d0/3.d0*zu_in_V(ix-3)+zu_in_V(ix+4)-zu_in_V(ix-4))/(-280.d0*dx)
    !dif2=(zu_in_V(ix+1)-2.d0*zu_in_V(ix)+zu_in_V(ix-1))/(dx*dx)
    !dif2=(-zu_in_V(ix+2)+16.d0*zu_in_V(ix+1)+16.d0*zu_in_V(ix-1)-zu_in_V(ix-2)-30.d0*zu_in_V(ix))/(12.d0*dx*dx)
    !dif2 = (-896.d0*zu_in_V(ix+1)-896.d0*zu_in_V(ix-1)+112.d0*zu_in_V(ix+2)+112.d0*zu_in_V(ix-2) &
    !      & -128.d0/9.d0*zu_in_V(ix+3)-128.d0/9.d0*zu_in_V(ix-3)+zu_in_V(ix+4)+zu_in_V(ix-4)+2.d0*7175.d0/9.d0*zu_in_V(ix))/(-560.d0*dx*dx)
    !hzu_V(ix)=0.5d0*(-dif2 -2.d0*zI*(At_in+k_in)*dif +(At_in+k_in)*(At_in+k_in)*zu_in_V(ix))+pot(ix)*zu_in_V(ix)



    v=(lap(1)*(zu_in_V(IDX(1))+zu_in_V(IDX(-1))) &
    & +lap(2)*(zu_in_V(IDX(2))+zu_in_V(IDX(-2))) &
    & +lap(3)*(zu_in_V(IDX(3))+zu_in_V(IDX(-3))) &
    & +lap(4)*(zu_in_V(IDX(4))+zu_in_V(IDX(-4))))
    w=(nab(1)*(zu_in_V(IDX(1))-zu_in_V(IDX(-1))) &
    & +nab(2)*(zu_in_V(IDX(2))-zu_in_V(IDX(-2))) &
    & +nab(3)*(zu_in_V(IDX(3))-zu_in_V(IDX(-3))) &
    & +nab(4)*(zu_in_V(IDX(4))-zu_in_V(IDX(-4))))

    hzu_V(ix) = pot(ix)*zu_in_V(ix) + 0.5d0*((k_in+At_in)*(k_in+At_in)-lap(0))*zu_in_V(ix) - 0.5d0 * v - zI * (k_in+At_in) * w 
  end do

END SUBROUTINE
