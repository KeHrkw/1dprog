SUBROUTINE ground_state()
use CONSTANTS
use WAVE_FUNC
use COEFF, only : eps
  implicit none
  integer :: iter,ik,ib
  real(8),dimension(:) :: hav_base_old(Nb),cur(Nb)
  integer,dimension(:) :: iflag_CG_conv(1:Nk)
  !open(8,file="./log_base.log")
  hav_base_old(:)=0d0
  do iter = 1, Niter, 1
    write(*,*) "---------------------------------------"
    write(*,*) "iter=",iter

    call CG_method(iflag_CG_conv)
    !call CG_method_fix()

    cur(:)=0d0
    do ib = 1, Nb, 1
      do ik = 1, Nk, 1
        k_in=k(ik)
        call zh_Velocity_operation(k_in,u(:,ik,ib),hzu(:,0))
        call Current_Velocity_operation(k_in,u(:,ik,ib),czu(:,0))
        cur(ib)=cur(ib)+real(sum(conjg(u(:,ik,ib))*czu(:,0))*dx)
        eps(ik,ib)=real(sum(conjg(u(:,ik,ib))*hzu(:,0))*dx)
      end do
      hav_base(ib)=sum(eps(:,ib))/real(Nk)
      cur(ib)=cur(ib)/real(Nk)
    end do
    write(*,*) "Energy of ik=1"
    do ib=1,3
      write(*,'(i,f11.6)',advance='no') ib,eps(1,ib)
    end do
    write(*,*)
    write(*,*) "Current"
    do ib=1,3
      write(*,'(i,e11.4)',advance='no') ib,cur(ib)
    end do
    write(*,*)
    write(*,*) "Total energy"
    do ib=1,Nb
      write(*,*) ib,hav_base(ib),hav_base(ib)-hav_base_old(ib)
    end do
    write(*,*) "---------------------------------------"
    if(sum(iflag_CG_conv(:))==Nk) then
      write(*,*) "----grand state calculation converge----"
      exit
    end if
    hav_base_old(:)=hav_base(:)
    iflag_CG_conv=0
  end do
  !close(8)

end subroutine
SUBROUTINE CG_method(iflag_CG_conv)
  !$ use omp_lib
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  integer,dimension(:),intent(out) :: iflag_CG_conv(1:Nk)

  complex(8),dimension(:,:) :: gk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: pk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: hpk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: xk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: hxk(0:Nx-1,0:NUMBER_THREADS-1)

  real(8) :: gkgk,pkpk,beta,xkxk, xkHxk,pkHpk,R
  real(8) :: hav_old
  complex(8) :: alpha,A,xkHpk,xkpk,B,Delta,C
  integer :: iscf, ib, jb, ik
  integer :: thr_id
  !$ double precision st, en
  thr_id=0
  iflag_CG_conv=0

  !$ st = omp_get_wtime()
  !write(8,*) "CG method"
  !write(8,*) "ik,ib,iscf,(xkHxk/xkxk-hav_old)"
  !$omp parallel private(thr_id)
  !$ thr_id=omp_get_thread_num()
  !$omp do private(ik,ib,jb,xkxk,xkHxk,hav_old,R,iscf,beta,gkgk,xkpk,pkpk,xkHpk,pkHpk,A,B,C,Delta,alpha)
  do ik = 1, Nk, 1
    do ib = 1, Nb, 1

      !Gram-Schmidt process
      do jb = 1, ib-1, 1
        u(:,ik,ib)=u(:,ik,ib)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*u(:,ik,ib))*dx
      end do


      xkxk=sum(conjg(u(:,ik,ib))*u(:,ik,ib))*dx
      xk(:,thr_id)=u(:,ik,ib)/sqrt(xkxk)

      !CG-method
      call zh_Velocity_operation(k(ik),xk(:,thr_id),hxk(:,thr_id))
      xkHxk=sum(conjg(xk(:,thr_id))*hxk(:,thr_id))*dx
      hav_old=xkHxk
      R=xkHxk/xkxk


      do iscf = 1, Nscf, 1

        gk(:,thr_id)=2.d0*(hxk(:,thr_id)-R*xk(:,thr_id))/xkxk
        !Gram-Schmidt process
        do jb = 1, ib-1, 1
          gk(:,thr_id)=gk(:,thr_id)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*gk(:,thr_id))*dx
        end do
        beta=sum(abs(gk(:,thr_id))**2)*dx
        if(iscf==1)then
          pk(:,thr_id)=-gk(:,thr_id)
        else
          pk(:,thr_id)=-gk(:,thr_id)+beta/gkgk*pk(:,thr_id)
        end if
        gkgk=beta

        xkpk=sum(conjg(xk(:,thr_id))*pk(:,thr_id))*dx
        pkpk=sum(abs(pk(:,thr_id))**2)*dx

        call zh_Velocity_operation(k(ik),pk(:,thr_id),hpk(:,thr_id))

        pkHpk=sum(conjg(pk(:,thr_id))*hpk(:,thr_id))*dx
        xkHpk=sum(conjg(xk(:,thr_id))*hpk(:,thr_id))*dx

        A=pkHpk*xkpk-xkHpk*pkpk
        B=pkHpk*xkxk-xkHxk*pkpk
        C=xkHpk*xkxk-xkHxk*xkpk

        Delta=B**2-4.d0*A*C
        if ( real(Delta)<0.d0 ) then
          write(*,'(<3>i,<6>e,a)') ik,ib,iscf,A,B,C,"ERE!"
          !write(8,*)A,B,C
          !write(8,*)alpha,beta,Delta
          !write(8,*)xkHxk,xkHpk
          !write(8,*)pkpk,gkgk,xkpk
          exit
        else
          alpha=(-B+sqrt(real(Delta)))/(2.d0*A)
        end if

         xk(:,thr_id)= xk(:,thr_id)+alpha* pk(:,thr_id)
        hxk(:,thr_id)=hxk(:,thr_id)+alpha*hpk(:,thr_id)

        !write(8,'(<3>i,<6>e,i)') ik,ib,iscf,xkxk,xkHxk,(R-hav_old),(gkgk*xkxk/(R**2)),Delta,thr_id

        if(gkgk*xkxk/(R**2)<4.d-20) then
          !write(*,*) ib,ik,"Converge"
          if(ib==Nb) iflag_CG_conv(ik)=1
          exit
        end if

        xkHxk=sum(conjg(xk(:,thr_id))*hxk(:,thr_id))*dx
        xkxk =sum(abs(xk(:,thr_id))**2)*dx
        R=xkHxk/xkxk

        !if(abs(R-hav_old) < 4.d-15) then
        !  write(*,'(<3>i,<2>e)') ik,ib,iter,(gkgk*xkxk/(R**2)),abs(R-hav_old)
        !exit
        !end if


        hav_old=R

      end do
      u(:,ik,ib)=xk(:,thr_id)/sqrt(xkxk)
    end do
  end do
  !$omp end do
  !$omp end parallel

  !$ en = omp_get_wtime()
  !$ print *, "Elapsed time [sec]:", en-st
END SUBROUTINE
SUBROUTINE CG_method_fix()
  !$ use omp_lib
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  complex(8),dimension(:,:) :: gk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: pk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: hpk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: xk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: hxk(0:Nx-1,0:NUMBER_THREADS-1)
  complex(8),dimension(:,:) :: pko(0:Nx-1,0:NUMBER_THREADS-1)

  real(8) :: gkgk,xkHxk,pkHpk,uk,s,ev
  complex(8) :: xkHpk,xkpk,cx,cp,zs
  integer :: iscf, ib, ibt, ik
  real(8),parameter :: delta_cg=1.d-15
  real(8) :: hav_old
  integer :: thr_id
  !write(8,*) "CG method"
  !write(8,*) "ik,ib,iscf,xkHxk,abs(ev-xkHxk)"

  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ib,ibt,s,xkHxk,iscf,uk,gkgk,xkHpk,pkHpk,ev,cx,cp,zs)
  do ik=1,Nk
    k_in=k(ik)
  do ib=1,Nb
    do ibt=1,ib-1
      s=sum(conjg(u(:,ik,ibt))*u(:,ik,ib))*dx
      u(0:Nx-1,ik,ib)=u(0:Nx-1,ik,ib)-u(0:Nx-1,ik,ibt)*s
    end do
    s=1.0d0/sqrt(sum(abs(u(:,ik,ib))**2)*dx)
    xk(0:Nx-1,thr_id)=u(0:Nx-1,ik,ib)*s

    call zh_Velocity_operation(k_in,xk(:,thr_id),hxk(:,thr_id))

    xkHxk=sum(conjg(xk(:,thr_id))*hxk(:,thr_id))*dx

    do iscf=1,Nscf
      gk(0:Nx-1,thr_id)=(hxk(0:Nx-1,thr_id)-xkHxk*xk(0:Nx-1,thr_id))
      do ibt=1,ib-1
        zs=sum(conjg(u(:,ik,ibt))*gk(:,thr_id))*dx
        gk(0:Nx-1,thr_id)=gk(0:Nx-1,thr_id)-u(0:Nx-1,ik,ibt)*zs
      end do
      s=sum(abs(gk(:,thr_id))**2)*dx

      select case (iscf)
      case(1)
        pk(0:Nx-1,thr_id)=gk(0:Nx-1,thr_id)
      case default
        uk=s/gkgk
        pk(0:Nx-1,thr_id)=gk(0:Nx-1,thr_id)+uk*pk(0:Nx-1,thr_id)
      end select
      gkgk=s

      zs=sum(conjg(xk(:,thr_id))*pk(:,thr_id))*dx
      pko(0:Nx-1,thr_id)=pk(0:Nx-1,thr_id)-xk(0:Nx-1,thr_id)*zs
      s=1.0d0/sqrt(sum(abs(pko(:,thr_id))**2)*dx)
      pko(0:Nx-1,thr_id)=pko(0:Nx-1,thr_id)*s

      call zh_Velocity_operation(k_in,pko(:,thr_id),hpk(:,thr_id))

      xkHpk=sum(conjg( xk(:,thr_id))*hpk(:,thr_id))*dx
      pkHpk=sum(conjg(pko(:,thr_id))*hpk(:,thr_id))*dx
      ev=0.5d0*((xkHxk+pkHpk)-sqrt((xkHxk-pkHpk)**2+4*abs(xkHpk)**2))
      cx=xkHpk/(ev-xkHxk)
      cp=1.d0/sqrt(1.d0+abs(cx)**2)
      cx=cx*cp

      !write(8,'(<3>i,<2>e)') ik,ib,iscf,xkHxk,abs(ev-xkHxk)
      if(abs(ev-xkHxk)<delta_cg) exit
       xk(0:Nx-1,thr_id)=cx* xk(0:Nx-1,thr_id)+cp*pko(0:Nx-1,thr_id)
      hxk(0:Nx-1,thr_id)=cx*hxk(0:Nx-1,thr_id)+cp*hpk(0:Nx-1,thr_id)
      xkHxk=sum(conjg(xk(:,thr_id))*hxk(:,thr_id))*dx
    enddo

    s=1.0d0/sqrt(sum(abs(xk(:,thr_id))**2)*dx)
    u(0:Nx-1,ik,ib)=xk(0:Nx-1,thr_id)*s
  enddo
  enddo
!$omp end parallel
END SUBROUTINE
SUBROUTINE simpson(x,y,zs)
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  complex(kind(0d0)),dimension(0:Nx-1),intent(in) :: x,y
  complex(kind(0d0)),intent(out) :: zs
  complex(kind(0d0)),dimension(0:Nx/4-1) :: A,B,C,D
  integer :: ix,m

  do ix = 0, Nx-1, 1
    m=mod(ix,4)
    select case(m)
    case(0)
      A(ix/4)=conjg(x(ix))*y(ix)
    case(1)
      B(ix/4)=conjg(x(ix))*y(ix)
    case(2)
      C(ix/4)=conjg(x(ix))*y(ix)
    case(3)
      D(ix/4)=conjg(x(ix))*y(ix)
    end select
  end do

  zs=sum(14.d0*A(:)+32.d0*B(:)+12.d0*C(:)+32.d0*D(:))*2.d0*dx/45.d0

END SUBROUTINE

Subroutine Gram_Schmidt_ompk()
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  integer :: ik,ib,ibt
  real(8) :: s
  complex(8) :: zov

  do ik=1,Nk
  do ib=1,Nb
    do ibt=1,ib-1
      zov=sum(conjg(u(:,ik,ibt))*u(:,ik,ib))*dx
      u(:,ik,ib)=u(:,ik,ib)-u(:,ik,ibt)*zov
    enddo
    s=sum(abs(u(:,ik,ib))**2)*dx
    u(:,ik,ib)=u(:,ik,ib)/sqrt(s)
  enddo
  enddo

End Subroutine
