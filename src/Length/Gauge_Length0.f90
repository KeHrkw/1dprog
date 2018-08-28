SUBROUTINE Gauge_Length()
  !Time propagation in Length gauge !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  !$ use omp_lib
  implicit none
  !$omp threadprivate(zu_in_L,hzu_L,czu_L)
  real(8),dimension(-LNk:RNk,Nbt) :: norm_l, hav_l, cur_l
  complex(kind(0d0)) :: zcoef, zuix
  complex(kind(0d0)), dimension(:,:,:) :: czu_out(0:Nx-1,-LNk:RNk,Nbt)
  integer :: iexp, ik, ib, ix

  norm_l(:,:)=0.d0
  cur_l(:,:)=0.d0
  hav_l(:,:)=0.d0

  !$omp parallel default(shared), private(ib,zcoef,iexp,ik)
  !$omp do
  do ib = 1, Nbt, 1
  !$ nthreads=omp_get_num_threads()

    select case(it)
    case(0)
    case default
      zcoef=1.d0
      zu_in_L(0:Nx-1,-LNk:RNk)=zu(0:Nx-1,:,ib)

      do iexp=1,Nexp,1
        zcoef=zcoef*(-zI)*dt/iexp
        call zh_Length_operation(ib)
        zu(:,:,ib)=zu(:,:,ib)+zcoef*hzu_L(:,:)

        zu_in_L(0:Nx-1,-LNk:RNk)=hzu_L(0:Nx-1,:)
      end do
    end select


    zu_in_L(0:Nx-1,-LNk:RNk)=zu(0:Nx-1,:,ib)

    call zh_Length_operation(ib)
    call Current_Length_operation()

    do ik = -LNk, RNk, 1
      czu_out(:,ik,ib)=czu_L(:,ik)
      norm_l(ik,ib) = abs(sum(conjg(zu(:,ik,ib))*zu(:,ik,ib))*dx)
      hav_l(ik,ib) = real(sum(conjg(zu(:,ik,ib))*hzu_L(:,ik))*dx)
      cur_l(ik,ib) = real(sum(conjg(zu(:,ik,ib))*czu_L(:,ik))*dx)
    end do
  end do
  !$omp end do
  !$omp end parallel
  if(mod(it,5000)==0)then
    !do ik = -LNk, RNk, 1
      !write(9,'(<2>i,<4>e)') it,ik,cur_l(ik,1),hav_l(ik,1),cur_l(ik,2),hav_l(ik,2)
    !end do
     do ix = 0, Nx-1, 1
        write(9,'(<2>i,<7>e)') it,ix,sum(abs(zu(ix,:,1))**2)/real(Nk),sum(conjg(zu(ix,:,1))*czu_out(ix,:,1))/real(Nk), &
                              & Et_in*ix*dx,sum(abs(zu(ix,:,2))**2)/real(Nk),sum(conjg(zu(ix,:,2))*czu_out(ix,:,1))/real(Nk)
     end do
    write(9,*)''
    write(9,*)''
  end if

  do ib = 1, Nbt, 1
    norm(ib)=sum(norm_l(:,ib))
    hav(ib,it)=sum(hav_l(:,ib))/real(Nk)
    cur(ib,it)=sum(cur_l(:,ib))/real(Nk)
  end do
END SUBROUTINE
