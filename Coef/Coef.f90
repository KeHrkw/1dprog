SUBROUTINE Coef_init()
  use CONSTANTS
  use WAVE_FUNC
  use COEFF
  use FD_K
  implicit none
  !complex(kind(0d0)) :: dif_x, dif_k
  complex(kind(0d0)) :: w,x
  integer :: ib,jb,ix,ik
# define IDX(dt) modx(ix+(dt)+Nx)

  write(*,*) "Norb : ",Norb
  udku(:,:,:)=0.d0
  udxu(:,:,:)=0.d0
  do ib = 1, Norb, 1
    zu_in_L(0:Nx-1,-LNk:RNk)=u(0:Nx-1,:,ib)
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
    !zu_in_L(-1,:)  =zu_in_L(Nx-1,:)
    !zu_in_L(Nx,:)  =zu_in_L(0,:)
    !zu_in_L(-2,:)  =zu_in_L(Nx-2,:)
    !zu_in_L(Nx+1,:)=zu_in_L(1,:)
    !zu_in_L(-3,:)  =zu_in_L(Nx-3,:)
    !zu_in_L(Nx+2,:)=zu_in_L(2,:)
    !zu_in_L(-4,:)  =zu_in_L(Nx-4,:)
    !zu_in_L(Nx+3,:)=zu_in_L(3,:)

    do jb = 1, Norb, 1

      do ik = -LNk, RNk, 1
        do ix = 0, Nx-1, 1
          !idf_x=(zu_in_L(ix+1,ik)-zu_in_L(ix-1,ik))/(2.d0*dx)
          !dif_x=(-zu_in_L(ix+2,ik) +8.d0*zu_in_L(ix+1,ik) -8.d0*zu_in_L(ix-1,ik)+zu_in_L(ix-2,ik))/(12.d0*dx)
          !dif_x = (-224.d0*zu_in_L(ix+1,ik)+224.d0*zu_in_L(ix-1,ik)+56.d0*zu_in_L(ix+2,ik)-56.d0*zu_in_L(ix-2,ik) &
          !      & -32.d0/3.d0*zu_in_L(ix+3,ik)+32.d0/3.d0*zu_in_L(ix-3,ik)+zu_in_L(ix+4,ik)-zu_in_L(ix-4,ik))/(-280.d0*dx)
          !dif_k =(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1))/(2.d0*dk)
          !dif_k =(-zu_in_L(ix,ik+2) +8.d0*zu_in_L(ix,ik+1) -8.d0*zu_in_L(ix,ik-1)+zu_in_L(ix,ik-2))/(12.d0*dk)
          !dif_k = (zu_in_L(ix,ik+3)-9.d0*zu_in_L(ix,ik+2)+45.d0*zu_in_L(ix,ik+1)-45.d0*zu_in_L(ix,ik-1)+9.d0*zu_in_L(ix,ik-2)-zu_in_L(ix,ik-3))/(60.d0*dk)
          !dif_k = (-224.d0*zu_in_L(ix,ik+1)+224.d0*zu_in_L(ix,ik-1)+56.d0*zu_in_L(ix,ik+2)-56.d0*zu_in_L(ix,ik-2) &
          !      & -32.d0/3.d0*zu_in_L(ix,ik+3)+32.d0/3.d0*zu_in_L(ix,ik-3)+zu_in_L(ix,ik+4)-zu_in_L(ix,ik-4))/(-280.d0*dk)
          !dif_k = (zu_in_L(ix,ik+5)-zu_in_L(ix,ik-5)-25.d0/2.d0*zu_in_L(ix,ik+4)+25.d0/2.d0*zu_in_L(ix,ik-4)+75.d0*zu_in_L(ix,ik+3)-75.d0*zu_in_L(ix,ik-3) &
          !      & -300.d0*zu_in_L(ix,ik+2)+300.d0*zu_in_L(ix,ik-2)+1050.d0*zu_in_L(ix,ik+1)-1050.d0*zu_in_L(ix,ik-1))/(2.d0*630.d0*dk)

          w=(nab(1)*(zu_in_L(IDX(1),ik)-zu_in_L(IDX(-1),ik)) &
          & +nab(2)*(zu_in_L(IDX(2),ik)-zu_in_L(IDX(-2),ik)) &
          & +nab(3)*(zu_in_L(IDX(3),ik)-zu_in_L(IDX(-3),ik)) &
          & +nab(4)*(zu_in_L(IDX(4),ik)-zu_in_L(IDX(-4),ik)))

          x=(nab_k(1)*(zu_in_L(ix,ik+1)-zu_in_L(ix,ik-1)) &
          & +nab_k(2)*(zu_in_L(ix,ik+2)-zu_in_L(ix,ik-2)) &
          & +nab_k(3)*(zu_in_L(ix,ik+3)-zu_in_L(ix,ik-3)) &
          & +nab_k(4)*(zu_in_L(ix,ik+4)-zu_in_L(ix,ik-4)) &
          & +nab_k(5)*(zu_in_L(ix,ik+5)-zu_in_L(ix,ik-5)))

          udku(ik,jb,ib)=udku(ik,jb,ib)+conjg(u(ix,ik,jb))*x*dx
          udxu(ik,jb,ib)=udxu(ik,jb,ib)+conjg(u(ix,ik,jb))*w*dx
        end do
      end do
    end do
  end do
  Cnm(:,:,:)=0.d0
  do ib = 1, Nbt, 1
    do jb = 1, Norb, 1
      if(ib==jb) Cnm(:,ib,jb)=1.d0
    end do
  end do
END SUBROUTINE
SUBROUTINE Coef_Length()
  use CONSTANTS
  use TD_CALC
  use COEFF
  use WAVE_FUNC
  implicit none
  real(8),dimension(-LNk:RNk,Nbt) :: norm_c, hav_c, cur_c
  complex(kind(0d0)) :: dum(-LNk:RNk,Nbt,Norb)
  complex(kind(0d0)) :: zcoef
  integer :: ib,jb,iexp,ik
  norm_c(:,:)=0.d0
  cur_c(:,:)=0.d0
  hav_c(:,:)=0.d0

  do ib = 1, Nbt, 1
    select case(it)
    case(0)
    case default
      zcoef=1.d0
      do jb = 1, Norb, 1
        Cnm_in(-LNk:RNk,jb)=Cnm(:,ib,jb)
      end do

      do iexp=1,Nexp,1
        zcoef=zcoef*(-zI)*dt/iexp

        call zh_Coef_Length()
        do jb = 1, Norb, 1
          Cnm(:,ib,jb)=Cnm(:,ib,jb)+zcoef*HCnm(:,jb)
          Cnm_in(-LNk:RNk,jb)=HCnm(-LNk:RNk,jb)
        end do
      end do
    end select

    do jb = 1, Norb, 1
      Cnm_in(-LNk:RNk,jb)=Cnm(:,ib,jb)
    end do
    call zh_Coef_Length()

    do ik = -LNk, RNk, 1
      norm_c(ik,ib) = abs(sum(conjg(Cnm(ik,ib,:))*Cnm(ik,ib,:)))
      hav_c(ik,ib) = real(sum(conjg(Cnm(ik,ib,:))*HCnm(ik,:)))
      do jb = 1, Norb, 1
        dum(ik,ib,jb) = sum(Cnm(ik,ib,:)*udxu(ik,jb,:))
      end do
      cur_c(ik,ib) = real(-zI*sum(conjg(Cnm(ik,ib,:))*dum(ik,ib,:)))+norm_c(ik,ib)*k(ik)
    end do
  end do
  do ib = 1, Nbt, 1
    norm(ib)=sum(norm_c(:,ib))
    hav(ib,it)=sum(hav_c(:,ib))/real(Nk)
    cur(ib,it)=sum(cur_c(:,ib))/real(Nk)
  end do

  if(mod(it,5000)==0)then
    do ib = 1, 1, 1
      do ik = -LNk, RNk, 1
        do jb = 1, Norb, 1
          write(11,'(<3>i,<2>e)') ib,ik,jb,Cnm(ik,ib,jb)
        end do
        write(11,*)""
      end do
      write(11,*)""
      write(11,*)""
    end do
  end if

END SUBROUTINE
SUBROUTINE zh_Coef_Length()
  use CONSTANTS
  use TD_CALC
  use COEFF
  use FD_K
  implicit none
  !complex(kind(0d0)) :: dif_k
  complex(kind(0d0)) :: x
  integer :: ik, jb
  Cnm_in(-LNk-1,:)=Cnm_in(RNk,:)
  Cnm_in(-LNk-2,:)=Cnm_in(RNk-1,:)
  Cnm_in(-LNk-3,:)=Cnm_in(RNk-2,:)
  Cnm_in(-LNk-4,:)=Cnm_in(RNk-3,:)
  Cnm_in(-LNk-5,:)=Cnm_in(RNk-4,:)
  Cnm_in(RNk+1,:)=Cnm_in(-LNk,:)
  Cnm_in(RNk+2,:)=Cnm_in(-LNk+1,:)
  Cnm_in(RNk+3,:)=Cnm_in(-LNk+2,:)
  Cnm_in(RNk+4,:)=Cnm_in(-LNk+3,:)
  Cnm_in(RNk+5,:)=Cnm_in(-LNk+4,:)
  do ik = -LNk, RNk, 1
    do jb = 1, Norb, 1
      !dif_k = (Cnm_in(ik+5,jb)-Cnm_in(ik-5,jb)-25.d0/2.d0*Cnm_in(ik+4,jb)+25.d0/2.d0*Cnm_in(ik-4,jb)+75.d0*Cnm_in(ik+3,jb)-75.d0*Cnm_in(ik-3,jb) &
      !      & -300.d0*Cnm_in(ik+2,jb)+300.d0*Cnm_in(ik-2,jb)+1050.d0*Cnm_in(ik+1,jb)-1050.d0*Cnm_in(ik-1,jb))/(2.d0*630.d0*dk)

        x=(nab_k(1)*(Cnm_in(ik+1,jb)-Cnm_in(ik-1,jb)) &
        & +nab_k(2)*(Cnm_in(ik+2,jb)-Cnm_in(ik-2,jb)) &
        & +nab_k(3)*(Cnm_in(ik+3,jb)-Cnm_in(ik-3,jb)) &
        & +nab_k(4)*(Cnm_in(ik+4,jb)-Cnm_in(ik-4,jb)) &
        & +nab_k(5)*(Cnm_in(ik+5,jb)-Cnm_in(ik-5,jb)))

      HCnm(ik,jb)=Cnm_in(ik,jb)*eps(ik,jb)+zI*Et_in*x+zI*Et_in*sum(Cnm_in(ik,:)*udku(ik,jb,:))
    end do
  end do
END SUBROUTINE
