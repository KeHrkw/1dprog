SUBROUTINE zh_Length_operation(ib)
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none
  !$omp threadprivate(zu_in_L,hzu_L,phase_in)
  integer,intent(in) :: ib
  integer :: ix, ik
  complex(kind(0d0)) :: dif_x, dif2, dif_k

  do ix = 0, Nx-1, 1
    zu_in_L(ix,-LNk-1)=phase(ix,ib,0)*zu_in_L(ix,RNk)
    zu_in_L(ix,RNk+1) =phase(ix,ib,1)*zu_in_L(ix,-LNk)
    zu_in_L(ix,-LNk-2)=phase(ix,ib,0)*zu_in_L(ix,RNk-1)
    zu_in_L(ix,RNk+2) =phase(ix,ib,1)*zu_in_L(ix,-LNk+1)
    zu_in_L(ix,-LNk-3)=phase(ix,ib,0)*zu_in_L(ix,RNk-2)
    zu_in_L(ix,RNk+3) =phase(ix,ib,1)*zu_in_L(ix,-LNk+2)
    zu_in_L(ix,-LNk-4)=phase(ix,ib,0)*zu_in_L(ix,RNk-3)
    zu_in_L(ix,RNk+4) =phase(ix,ib,1)*zu_in_L(ix,-LNk+3)
    zu_in_L(ix,-LNk-5)=phase(ix,ib,0)*zu_in_L(ix,RNk-4)
    zu_in_L(ix,RNk+5) =phase(ix,ib,1)*zu_in_L(ix,-LNk+4)
  end do
  zu_in_L(-1,:)  =zu_in_L(Nx-1,:)
  zu_in_L(Nx,:)  =zu_in_L(0,:)
  zu_in_L(-2,:)  =zu_in_L(Nx-2,:)
  zu_in_L(Nx+1,:)=zu_in_L(1,:)
  zu_in_L(-3,:)  =zu_in_L(Nx-3,:)
  zu_in_L(Nx+2,:)=zu_in_L(2,:)
  zu_in_L(-4,:)  =zu_in_L(Nx-4,:)
  zu_in_L(Nx+3,:)=zu_in_L(3,:)
  !$omp parallel private(dif2, dif_x, dif_k)
  !$omp do
  do ik = -LNk, RNk, 1
    do ix = 0, Nx-1, 1
      !idf_x=(zu_in_L(ix+1,ik)-zu_in_L(ix-1,ik))/(2.d0*dx)
      !dif_x=(-zu_in_L(ix+2,ik) +8.d0*zu_in_L(ix+1,ik) -8.d0*zu_in_L(ix-1,ik)+zu_in_L(ix-2,ik))/(12.d0*dx)
      dif_x = (-224.d0*zu_in_L(ix+1,ik)+224.d0*zu_in_L(ix-1,ik)+56.d0*zu_in_L(ix+2,ik)-56.d0*zu_in_L(ix-2,ik) &
            & -32.d0/3.d0*zu_in_L(ix+3,ik)+32.d0/3.d0*zu_in_L(ix-3,ik)+zu_in_L(ix+4,ik)-zu_in_L(ix-4,ik))/(-280.d0*dx)
      !dif2 = (zu_in_L(ix+1,ik)-2.d0*zu_in_L(ix,ik)+zu_in_L(ix-1,ik))/(dx*dx)
      !dif2 = (-zu_in_L(ix+2,ik)+16.d0*zu_in_L(ix+1,ik)+16.d0*zu_in_L(ix-1,ik)-zu_in_L(ix-2,ik)-30.d0*zu_in_L(ix,ik))/(12.d0*dx*dx)
      dif2 = (-896.d0*zu_in_L(ix+1,ik)-896.d0*zu_in_L(ix-1,ik)+112.d0*zu_in_L(ix+2,ik)+112.d0*zu_in_L(ix-2,ik) &
            & -128.d0/9.d0*zu_in_L(ix+3,ik)-128.d0/9.d0*zu_in_L(ix-3,ik)+zu_in_L(ix+4,ik)+zu_in_L(ix-4,ik)+2.d0*7175.d0/9.d0*zu_in_L(ix,ik))/(-560.d0*dx*dx)

      !dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk)
      !dif_k =(-zu_in_L(ix,ik+2) +8.d0*zu_in_L(ix,ik+1) -8.d0*zu_in_L(ix,ik-1)+zu_in_L(ix,ik-2))/(12.d0*dk)
      !dif_k = (zu_in_L(ix,ik+3)-9.d0*zu_in_L(ix,ik+2)+45.d0*zu_in_L(ix,ik+1)-45.d0*zu_in_L(ix,ik-1)+9.d0*zu_in_L(ix,ik-2)-zu_in_L(ix,ik-3))/(60.d0*dk)
      !dif_k = (-224.d0*zu_in_L(ix,ik+1)+224.d0*zu_in_L(ix,ik-1)+56.d0*zu_in_L(ix,ik+2)-56.d0*zu_in_L(ix,ik-2) &
      !      & -32.d0/3.d0*zu_in_L(ix,ik+3)+32.d0/3.d0*zu_in_L(ix,ik-3)+zu_in_L(ix,ik+4)-zu_in_L(ix,ik-4))/(-280.d0*dk)
      dif_k = (zu_in_L(ix,ik+5)-zu_in_L(ix,ik-5)-25.d0/2.d0*zu_in_L(ix,ik+4)+25.d0/2.d0*zu_in_L(ix,ik-4)+75.d0*zu_in_L(ix,ik+3)-75.d0*zu_in_L(ix,ik-3) &
            & -300.d0*zu_in_L(ix,ik+2)+300.d0*zu_in_L(ix,ik-2)+1050.d0*zu_in_L(ix,ik+1)-1050.d0*zu_in_L(ix,ik-1))/(2.d0*630.d0*dk)

      !if(ik < LNk_1 .AND. ik > -LNk_1) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk_1)
      !else if ( (ik > LNk_1) .AND. (ik < LNk_1+LNk_2) ) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk_2)
      !else if ( (ik < -LNk_1) .AND. (ik > -LNk_1-LNk_2)) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk_2)
      !else if ( (ik > LNk_1+LNk_2 )) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk_1)
      !else if ( (ik < -LNk_1-LNk_2 )) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk_1)
      !else if ( (ik==-LNk_1-LNk_2) .OR. (ik==LNk_1+LNk_2) .OR. (ik==-LNk_1) .OR. (ik==LNk_1)  ) then
      !  dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(dk_2+dk_1)
      !end if

      hzu_L(ix,ik)=0.5d0*(-dif2 -2.d0*zI*(k(ik))*dif_x +(k(ik))*(k(ik))*zu_in_L(ix,ik)) +pot(ix)*zu_in_L(ix,ik)+zI*Et_in*dif_k
    end do
  end do
  !$omp end do
  !$omp end parallel
END SUBROUTINE
