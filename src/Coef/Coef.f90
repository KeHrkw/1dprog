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
    zu_in_L(0:Nx-1,1:Nk)=u(0:Nx-1,:,ib)
    do ix = 0, Nx-1, 1
      zu_in_L(ix,0)=phase(ix,ib,0)*zu_in_L(ix,Nk)
      zu_in_L(ix,-1)=phase(ix,ib,0)*zu_in_L(ix,Nk-1)
      zu_in_L(ix,-2)=phase(ix,ib,0)*zu_in_L(ix,Nk-2)
      zu_in_L(ix,-3)=phase(ix,ib,0)*zu_in_L(ix,Nk-3)
      zu_in_L(ix,-4)=phase(ix,ib,0)*zu_in_L(ix,Nk-4)
      zu_in_L(ix,Nk+1) =phase(ix,ib,1)*zu_in_L(ix,1)
      zu_in_L(ix,Nk+2) =phase(ix,ib,1)*zu_in_L(ix,2)
      zu_in_L(ix,Nk+3) =phase(ix,ib,1)*zu_in_L(ix,3)
      zu_in_L(ix,Nk+4) =phase(ix,ib,1)*zu_in_L(ix,4)
      zu_in_L(ix,Nk+5) =phase(ix,ib,1)*zu_in_L(ix,5)
    end do

    do jb = 1, Norb, 1

      do ik = 1, Nk, 1
        do ix = 0, Nx-1, 1
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
SUBROUTINE Coef_Length(it,Et_in)
  use CONSTANTS
  use TD_CALC
  use COEFF
  use WAVE_FUNC
  implicit none
  integer,intent(in) :: it
  real(8),intent(in) :: Et_in
  real(8),dimension(1:Nk,Nbt) :: norm_c, hav_c, cur_c
  complex(kind(0d0)) :: dum(1:Nk,Nbt,Norb)
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
        Cnm_in(1:Nk,jb)=Cnm(:,ib,jb)
      end do

      do iexp=1,Nexp,1
        zcoef=zcoef*(-zI)*dt/iexp

        call zh_Coef_Length(Et_in)
        do jb = 1, Norb, 1
          Cnm(:,ib,jb)=Cnm(:,ib,jb)+zcoef*HCnm(:,jb)
          Cnm_in(1:Nk,jb)=HCnm(1:Nk,jb)
        end do
      end do
    end select

    do jb = 1, Norb, 1
      Cnm_in(1:Nk,jb)=Cnm(:,ib,jb)
    end do
    call zh_Coef_Length(Et_in)

    do ik = 1, Nk, 1
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
      do ik = 1, Nk, 1
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
SUBROUTINE zh_Coef_Length(Et_in)
  use CONSTANTS
  use TD_CALC
  use COEFF
  use FD_K
  implicit none
  !complex(kind(0d0)) :: dif_k
  real(8),intent(in) :: Et_in
  complex(kind(0d0)) :: x
  integer :: ik, jb
  Cnm_in(0,:)=Cnm_in(Nk,:)
  Cnm_in(-1,:)=Cnm_in(Nk-1,:)
  Cnm_in(-2,:)=Cnm_in(Nk-2,:)
  Cnm_in(-3,:)=Cnm_in(Nk-3,:)
  Cnm_in(-4,:)=Cnm_in(Nk-4,:)
  Cnm_in(Nk+1,:)=Cnm_in(1,:)
  Cnm_in(Nk+2,:)=Cnm_in(2,:)
  Cnm_in(Nk+3,:)=Cnm_in(3,:)
  Cnm_in(Nk+4,:)=Cnm_in(4,:)
  Cnm_in(Nk+5,:)=Cnm_in(5,:)
  do ik = 1, Nk, 1
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
