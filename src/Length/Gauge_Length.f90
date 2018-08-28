SUBROUTINE Gauge_Length()
  !Time propagation in velocity gauge!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  use FD_K
  !$ use omp_lib
  implicit none
  real(8),dimension(-LNk:RNk,Nbt) :: norm_v, hav_v, cur_v
  complex(kind(0d0)) :: zcoef
  complex(kind(0d0)), dimension(:,:,:) :: czu_out(0:Nx-1,-LNk:RNk,Nbt)
  integer, save :: iexp, ib
  integer :: ik, ix
  !$omp threadprivate(czu_V,hzu_V,zu_in_V,k_in,dk_zu,zu_in_L)

  norm_v(:,:)=0.d0
  cur_v(:,:)=0.d0
  hav_v(:,:)=0.d0

  !$omp parallel default(shared), private(ib,zcoef,iexp,ik)
  !$omp do
  do ib=1, Nbt, 1
      !$ nthreads=omp_get_num_threads()
    select case(it)
    case(0)
    case default
      zu_in_L(0:Nx-1,-LNk:RNk)=zu(0:Nx-1,-LNk:RNk,ib)
      zcoef=1.d0

      do iexp=1,Nexp,1
        zcoef=zcoef*(-zI)*dt/iexp

        call dev_k(ib)

        do ik = -LNk,RNk
          k_in=k(ik)
          zu_in_V(0:Nx-1)=zu_in_L(0:Nx-1,ik)
          call zh_Velocity_operation()

          hzu_V(:)=hzu_V(:)+zI*Et_in*dk_zu(:,ik)

          zu(:,ik,ib)=zu(:,ik,ib)+zcoef*hzu_V(:)
          zu_in_L(:,ik)=hzu_v(:)
        end do
      end do
    end select

    zu_in_L(0:Nx-1,-LNk:RNk)=zu(0:Nx-1,-LNk:RNk,ib)
    call dev_k(ib)

    do ik=-LNk,RNk
      k_in=k(ik)
      zu_in_V(0:Nx-1)=zu(0:Nx-1,ik,ib)
      call zh_Velocity_operation()
      hzu_V(:)=hzu_V(:)+zI*Et_in*dk_zu(:,ik)

      call Current_Velocity_operation()

      czu_out(:,ik,ib)=czu_V(:)
      norm_v(ik,ib) = abs(sum(conjg(zu(:,ik,ib))*zu(:,ik,ib))*dx)
      hav_v(ik,ib) = real(sum(conjg(zu(:,ik,ib))*hzu_V(:))*dx)
      cur_v(ik,ib) = real(sum(conjg(zu(:,ik,ib))*czu_V(:))*dx)
    end do
  end do
  !$omp end do
  !$omp end parallel

  if(mod(it,5000)==0)then
    !do ik = -LNk+1, LNk, 1
    !  write(9,'(<2>i,<2>e)') it,ik,cur_v(ik,1),hav_v(ik,1)
    !end do
    do ix = 0, Nx-1, 1
       write(9,'(<2>i,<7>e)') it,ix,sum(abs(zu(ix,:,1))**2)/real(Nk),sum(conjg(zu(ix,:,1))*czu_out(ix,:,1))/real(Nk), &
                             & Et(it)*ix*dx,sum(abs(zu(ix,:,2))**2)/real(Nk),sum(conjg(zu(ix,:,2))*czu_out(ix,:,1))/real(Nk)
    end do
    write(9,*)''
    write(9,*)''
  end if
  do ib = 1, Nbt, 1
    norm(ib)=sum(norm_v(:,ib))
    hav(ib,it)=sum(hav_v(:,ib))/real(Nk)
    cur(ib,it)=sum(cur_v(:,ib))/real(Nk)
    !norm(ib)=(sum(norm_v(-LNk_1:LNk_1,ib))+sum(norm_v(LNk_1+1:LNk_1+LNk_2*(LNk_part-2),ib))*real(LNk_1/LNk_2)+sum(norm_v(LNk_1+LNk_2*(LNk_part-2)+1:RNk,ib)) &
    !          & +sum(norm_v(-LNk_1-LNk_2*(LNk_part-2):-LNk_1-1,ib))*real(LNk_1/LNk_2)+sum(norm_v(-LNk:-LNk_1-LNk_2*(LNk_part-2)-1,ib)) )/real(LNk_1*LNk_part*2+1)
    !hav(ib,it)=(sum(hav_v(-LNk_1:LNk_1,ib))+sum(hav_v(LNk_1+1:LNk_1+LNk_2*(LNk_part-2),ib))*real(LNk_1/LNk_2)+sum(hav_v(LNk_1+LNk_2*(LNk_part-2)+1:RNk,ib)) &
    !          & +sum(hav_v(-LNk_1-LNk_2*(LNk_part-2):-LNk_1-1,ib))*real(LNk_1/LNk_2)+sum(hav_v(-LNk:-LNk_1-LNk_2*(LNk_part-2)-1,ib)) )/real(LNk_1*LNk_part*2+1)
    !cur(ib,it)=(sum(cur_v(-LNk_1:LNk_1,ib))+sum(cur_v(LNk_1+1:LNk_1+LNk_2*(LNk_part-2),ib))*real(LNk_1/LNk_2)+sum(cur_v(LNk_1+LNk_2*(LNk_part-2)+1:RNk,ib)) &
    !          & +sum(cur_v(-LNk_1-LNk_2*(LNk_part-2):-LNk_1-1,ib))*real(LNk_1/LNk_2)+sum(cur_v(-LNk:-LNk_1-LNk_2*(LNk_part-2)-1,ib)) )/real(LNk_1*LNk_part*2+1)
  end do

END SUBROUTINE
SUBROUTINE dev_k(ib)
  use CONSTANTS
  use FD_K
  !$ use omp_lib
  implicit none
  complex(8),dimension(0:Nx-1,-LNk-Nd_k:RNk+Nd_k) :: zu_L
  integer,intent(in) :: ib
  !complex(8),dimension(0:Nx-1,-LNk:RNk),intent(in) :: A
  !complex(8),dimension(0:Nx-1,-LNk:RNk),intent(out) :: B
  integer :: ik,ix
  do ix = 0, Nx-1, 1
    zu_L(ix,-LNk:RNk)=zu_in_L(ix,-LNk:RNk)
    zu_L(ix,-LNk-1)=phase(ix,ib,0)*zu_in_L(ix,RNk)
    zu_L(ix,-LNk-2)=phase(ix,ib,0)*zu_in_L(ix,RNk-1)
    zu_L(ix,-LNk-3)=phase(ix,ib,0)*zu_in_L(ix,RNk-2)
    zu_L(ix,-LNk-4)=phase(ix,ib,0)*zu_in_L(ix,RNk-3)
    zu_L(ix,-LNk-5)=phase(ix,ib,0)*zu_in_L(ix,RNk-4)
    zu_L(ix,RNk+1) =phase(ix,ib,1)*zu_in_L(ix,-LNk)
    zu_L(ix,RNk+2) =phase(ix,ib,1)*zu_in_L(ix,-LNk+1)
    zu_L(ix,RNk+3) =phase(ix,ib,1)*zu_in_L(ix,-LNk+2)
    zu_L(ix,RNk+4) =phase(ix,ib,1)*zu_in_L(ix,-LNk+3)
    zu_L(ix,RNk+5) =phase(ix,ib,1)*zu_in_L(ix,-LNk+4)
  end do

  do ik = -LNk, RNk, 1
    do ix = 0, Nx-1, 1
      dk_zu(ix,ik) = (nab_k(1)*(zu_L(ix,ik+1)-zu_L(ix,ik-1)) &
      &              + nab_k(2)*(zu_L(ix,ik+2)-zu_L(ix,ik-2)) &
      &              + nab_k(3)*(zu_L(ix,ik+3)-zu_L(ix,ik-3)) &
      &              + nab_k(4)*(zu_L(ix,ik+4)-zu_L(ix,ik-4)) &
      &              + nab_k(5)*(zu_L(ix,ik+5)-zu_L(ix,ik-5)))
    end do
  end do

  return
END SUBROUTINE dev_k
