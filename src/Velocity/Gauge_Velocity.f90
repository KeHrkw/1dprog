SUBROUTINE Gauge_Velocity(it,At_in)
  !Time propagation in velocity gauge!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  !$ use omp_lib
  implicit none
  integer,intent(in) :: it
  real(8),intent(in) :: At_in
  real(8),dimension(1:Nk,Nbt) :: norm_v, hav_v, cur_v
  complex(kind(0d0)) :: zcoef
  complex(kind(0d0)), dimension(:,:,:) :: czu_out(0:Nx-1,1:Nk,Nbt)
  integer, save :: iexp, ib
  integer :: ik, ix
  integer :: thr_id
  thr_id=0
  norm_v(:,:)=0.d0
  cur_v(:,:)=0.d0
  hav_v(:,:)=0.d0

  !$omp parallel default(shared), private(thr_id)
  !$ thr_id=omp_get_thread_num()
  !$omp do private(ik,ib,iexp,zcoef,k_in)
  do ik = 1, Nk, 1

    k_in=k(ik)+At_in
    do ib=1, Nbt, 1

      select case(it)
      case(0)
      case default
        zu_in(0:Nx-1,thr_id)=zu(0:Nx-1,ik,ib)
        zcoef=1.d0

        do iexp=1,Nexp,1
          zcoef=zcoef*(-zI)*dt/iexp
          call zh_Velocity_operation(k_in,zu_in(:,thr_id),hzu(:,thr_id))
          zu(:,ik,ib)=zu(:,ik,ib)+zcoef*hzu(:,thr_id)

          zu_in(0:Nx-1,thr_id)=hzu(0:Nx-1,thr_id)
        end do

      end select



      call zh_Velocity_operation(k_in,zu(:,ik,ib),hzu(:,thr_id))
      call Current_Velocity_operation(k_in,zu(:,ik,ib),czu(:,thr_id))

      czu_out(:,ik,ib)=czu(:,thr_id)
      norm_v(ik,ib) = abs(sum(conjg(zu(:,ik,ib))*zu(:,ik,ib))*dx)
      hav_v(ik,ib) = real(sum(conjg(zu(:,ik,ib))*hzu(:,thr_id))*dx)
      cur_v(ik,ib) = real(sum(conjg(zu(:,ik,ib))*czu(:,thr_id))*dx)
    end do
  end do
  !$omp end do
  !$omp end parallel

  do ib = 1, Nbt, 1
    norm(ib,it)=sum(norm_v(:,ib))
     hav(ib,it)=sum(hav_v(:,ib))/real(Nk)
     cur(ib,it)=sum(cur_v(:,ib))/real(Nk)
  end do

END SUBROUTINE
