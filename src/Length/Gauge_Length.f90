SUBROUTINE Gauge_Length(it,Et_in)
  !Time propagation in velocity gauge!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  use FD_K
  !$ use omp_lib
  implicit none
  integer,intent(in) :: it
  real(8),intent(in) :: Et_in
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
  do ib=1, Nbt, 1
    select case(it)
    case(0)
    case default
      zu_in_L(0:Nx-1,1:Nk,thr_id)=zu(0:Nx-1,1:Nk,ib)
      zcoef=1.d0

      do iexp=1,Nexp,1
        zcoef=zcoef*(-zI)*dt/iexp

        call dev_k(ib,zu_in_L(:,:,thr_id),dk_zu(:,:,thr_id))

        do ik = 1, Nk
          k_in=k(ik)
          call zh_Velocity_operation(k_in,zu_in_L(:,ik,thr_id),hzu(:,thr_id))
          zu_in_L(:,ik,thr_id)=hzu(:,thr_id)+zI*Et_in*dk_zu(:,ik,thr_id)
          zu(:,ik,ib)=zu(:,ik,ib)+zcoef*zu_in_L(:,ik,thr_id)
        end do
      end do
    end select

    call dev_k(ib,zu(:,:,ib),dk_zu(:,:,thr_id))

    do ik=1,Nk
      k_in=k(ik)
      call zh_Velocity_operation(k_in,zu(:,ik,ib),hzu(:,thr_id))
      hzu(:,thr_id)=hzu(:,thr_id)+zI*Et_in*dk_zu(:,ik,thr_id)

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
SUBROUTINE dev_k(ib,A,B)
  use CONSTANTS
  use FD_K
  implicit none
  complex(8),dimension(0:Nx-1,1:Nk),intent(in)  :: A
  complex(8),dimension(0:Nx-1,1:Nk),intent(out) :: B
  complex(8),dimension(0:Nx-1,1-Nd_k:Nk+Nd_k) :: zu_L
  integer,intent(in) :: ib
  integer :: ik,ix
  do ix = 0, Nx-1, 1
    zu_L(ix,1:Nk)=A(ix,1:Nk)
    zu_L(ix, 0)=phase(ix,ib,0)*A(ix,Nk)
    zu_L(ix,-1)=phase(ix,ib,0)*A(ix,Nk-1)
    zu_L(ix,-2)=phase(ix,ib,0)*A(ix,Nk-2)
    zu_L(ix,-3)=phase(ix,ib,0)*A(ix,Nk-3)
    zu_L(ix,-4)=phase(ix,ib,0)*A(ix,Nk-4)
    zu_L(ix,Nk+1) =phase(ix,ib,1)*A(ix,1)
    zu_L(ix,Nk+2) =phase(ix,ib,1)*A(ix,2)
    zu_L(ix,Nk+3) =phase(ix,ib,1)*A(ix,3)
    zu_L(ix,Nk+4) =phase(ix,ib,1)*A(ix,4)
    zu_L(ix,Nk+5) =phase(ix,ib,1)*A(ix,5)
  end do

  do ik = 1, Nk, 1
    do ix = 0, Nx-1, 1
      B(ix,ik) = (nab_k(1)*(zu_L(ix,ik+1)-zu_L(ix,ik-1)) &
      &              + nab_k(2)*(zu_L(ix,ik+2)-zu_L(ix,ik-2)) &
      &              + nab_k(3)*(zu_L(ix,ik+3)-zu_L(ix,ik-3)) &
      &              + nab_k(4)*(zu_L(ix,ik+4)-zu_L(ix,ik-4)) &
      &              + nab_k(5)*(zu_L(ix,ik+5)-zu_L(ix,ik-5)))
    end do
  end do

  return
END SUBROUTINE dev_k
