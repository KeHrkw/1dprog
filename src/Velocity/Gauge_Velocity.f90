SUBROUTINE Gauge_Velocity()
  !Time propagation in velocity gauge!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  !$ use omp_lib
  implicit none
  real(8),dimension(-LNk:RNk,Nbt) :: norm_v, hav_v, cur_v
  complex(kind(0d0)) :: zcoef
  complex(kind(0d0)), dimension(:,:,:) :: czu_out(0:Nx-1,-LNk:RNk,Nbt)
  integer, save :: iexp, ib
  integer :: ik, ix
  !$omp threadprivate(czu_V,hzu_V,zu_in_V,k_in,At_in)

  norm_v(:,:)=0.d0
  cur_v(:,:)=0.d0
  hav_v(:,:)=0.d0

  !$omp parallel default(shared), private(ik,ib,iexp,zcoef)
  !$omp do
  do ik = -LNk, RNk, 1
  !$ nthreads=omp_get_num_threads()
    k_in=k(ik)
    do ib=1, Nbt, 1

      select case(it)
      case(0)
      case default
        zu_in_V(0:Nx-1)=zu(0:Nx-1,ik,ib)
        zcoef=1.d0

        do iexp=1,Nexp,1
          zcoef=zcoef*(-zI)*dt/iexp
          call zh_Velocity_operation()
          zu(:,ik,ib)=zu(:,ik,ib)+zcoef*hzu_V(:)

          zu_in_V(0:Nx-1)=hzu_V(0:Nx-1)
        end do

      end select

      zu_in_V(0:Nx-1)=zu(0:Nx-1,ik,ib)


      call zh_Velocity_operation()
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
