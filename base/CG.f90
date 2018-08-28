SUBROUTINE ground_state()
use CONSTANTS
  implicit none
  integer :: iter
  open(8,file="./output/log_base.dat")

  do iter = 1, Niter, 1
    write(*,*) "---------------------------------"
    write(*,*) "iter=",iter
    call CG_method()
    write(*,*) "---------------------------------"
  end do
  close(8)

end subroutine
SUBROUTINE CG_method()
  use CONSTANTS
  use WAVE_FUNC
  use COEFF, only : eps
  implicit none
  complex(kind(0d0)),dimension(:) :: gk(0:Nx-1), pk(0:Nx-1), hpk(0:Nx-1), xk(0:Nx-1), hxk(0:Nx-1)

  real(8) :: gkgk,pkpk,beta,ev,xkxk, xkHxk,pkHpk,R,s
  real(8) :: hav_old
  complex(kind(0d0)) :: alpha,A,xkHpk,xkpk,B,Delta,C
  integer :: iscf, ib, jb, ik

  write(8,*) "CG method"
  write(8,*) "ik,ib,iscf,(xkHxk/xkxk-hav_old)"
  do ik = -LNk, RNk, 1
    k_in=k(ik)
    do ib = 1, Nb, 1

      !Gram-Schmidt process
      do jb = 1, ib-1, 1
        u(:,ik,ib)=u(:,ik,ib)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*u(:,ik,ib))*dx
      end do


      xkxk=sum(conjg(u(:,ik,ib))*u(:,ik,ib))*dx
      xk(:)=u(:,ik,ib)/sqrt(xkxk)

      !CG-method
      u_in(0:Nx-1)=xk(0:Nx-1)
      call h0_operation()
      hxk(:)=hu(:)
      xkHxk=sum(conjg(xk(:))*hxk(:))*dx
      hav_old=xkHxk
      R=xkHxk/xkxk


      do iscf = 1, Nscf, 1

        gk(:)=2.d0*(hxk(:)-R*xk(:))/xkxk
        !Gram-Schmidt process
        do jb = 1, ib-1, 1
          gk(:)=gk(:)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*gk(:))*dx
        end do
        beta=sum(abs(gk(:))**2)*dx
        if(iscf==1)then
          pk(:)=-gk(:)
        else
          pk(:)=-gk(:)+beta/gkgk*pk(:)
        end if
        gkgk=beta

        xkpk=sum(conjg(xk(:))*pk(:))*dx
        pkpk=sum(abs(pk(:))**2)*dx

        u_in(0:Nx-1)=pk(0:Nx-1)
        call h0_operation()
        hpk(:)=hu(:)

        pkHpk=sum(conjg(pk(:))*hpk(:))*dx
        xkHpk=sum(conjg(xk(:))*hpk(:))*dx

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

        xk(:)=xk(:)+alpha*pk(:)
        hxk(:)=hxk(:)+alpha*hpk(:)

        write(8,'(<3>i,<6>e)') ik,ib,iscf,xkxk,xkHxk,(R-hav_old),(gkgk*xkxk/(R**2)),Delta

        if(gkgk*xkxk/(R**2)<4.d-15) exit


        xkHxk=sum(conjg(xk(:))*hxk(:))*dx
        xkxk =sum(abs(xk(:))**2)*dx
        R=xkHxk/xkxk

        !if(abs(R-hav_old) < 4.d-15) then
        !  write(*,'(<3>i,<2>e)') ik,ib,iter,(gkgk*xkxk/(R**2)),abs(R-hav_old)
        !exit
        !end if


        hav_old=R

      end do
      u(:,ik,ib)=xk(:)/sqrt(xkxk)
    end do
  end do
  do ib = 1, Nb, 1
    do ik = -LNk, RNk, 1
      k_in=k(ik)
      u_in(0:Nx-1)=u(0:Nx-1,ik,ib)
      call h0_operation()
      eps(ik,ib)=real(sum(conjg(u(:,ik,ib))*hu(:))*dx)
    end do
    hav_base(ib)=sum(eps(:,ib))/real(Nk)
  end do
  write(*,*) "Energy of k=1"
  write(*,*) eps(1,:)
  write(*,*) "Total energy"
  write(*,*) hav_base(:)
END SUBROUTINE
SUBROUTINE CG_method_fix()
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  !$omp threadprivate(u_in,hu,k_in)
  complex(kind(0d0)),dimension(0:Nx-1) :: gk, pk, hpk, xk, hxk, pko

  real(8) :: gkgk,xkHxk,pkHpk,uk,s,ev
  complex(kind(0d0)) :: xkHpk,xkpk,cx,cp,zs
  integer :: iter, ib, jb, ik
  real(8),parameter :: delta_cg=1.d-7
  real(8) :: hav_old
  !write(8,*) "CG method"
  !write(8,*) "ik,ib,iter,abs(alpha-hav_old),(gkgk/(xkHxk**2))"
!$omp parallel
!$omp do private(ik,ib,jb,s,zs,xkHxk,iter,uk,gkgk,xkHpk,pkHpk,cx,cp,hav_old,ev,xkpk,gk, pk, hpk, xk, hxk, pko)
  do ik = -LNk, RNk, 1
    k_in=k(ik)
    do ib = 1, Nb, 1


      do jb=1,ib-1
        s=sum(conjg(u(:,ik,jb))*u(:,ik,ib))*dx
        u(0:Nx-1,ik,ib)=u(0:Nx-1,ik,ib)-u(0:Nx-1,ik,jb)*s
      end do
      s=1.0d0/sqrt(sum(abs(u(:,ik,ib))**2)*dx)
      xk(0:Nx-1)=u(0:Nx-1,ik,ib)*s

      u_in(0:Nx-1)=xk(0:Nx-1)
      call h0_operation()
      hxk(0:Nx-1)=hu(0:Nx-1)

      xkHxk=sum(conjg(xk(:))*hxk(:))*dx
      hav_old=xkHxk
      do iter=1,Niter

        gk(0:Nx-1)=(hxk(0:Nx-1)-xkHxk*xk(0:Nx-1))
        do jb=1,ib-1
          zs=sum(conjg(u(:,ik,jb))*gk(:))*dx
          gk(0:Nx-1)=gk(0:Nx-1)-u(0:Nx-1,ik,jb)*zs
        end do
        s=sum(abs(gk(:))**2)*dx

        select case (iter)
        case(1)
          pk(0:Nx-1)=gk(0:Nx-1)
        case default
          uk=s/gkgk
          pk(0:Nx-1)=gk(0:Nx-1)+uk*pk(0:Nx-1)
        end select
        gkgk=s

        zs=sum(conjg(xk(:))*pk(:))*dx
        pko(0:Nx-1)=pk(0:Nx-1)-xk(0:Nx-1)*zs
        s=1.0d0/sqrt(sum(abs(pko(:))**2)*dx)
        pko(0:Nx-1)=pko(0:Nx-1)*s

        u_in(0:Nx-1)=pko(0:Nx-1)
        call h0_operation()
        hpk(0:Nx-1)=hu(0:Nx-1)

        xkHpk=sum(conjg( xk(:))*hpk(:))*dx
        pkHpk=sum(conjg(pko(:))*hpk(:))*dx
        ev=0.5d0*((xkHxk+pkHpk)-sqrt((xkHxk-pkHpk)**2+4*abs(xkHpk)**2))
        cx=xkHpk/(ev-xkHxk)
        cp=1.d0/sqrt(1.d0+abs(cx)**2)
        cx=cx*cp
        if(abs(ev-xkHxk)<delta_cg) exit
        if(iter==Niter) exit
         xk(0:Nx-1)=cx* xk(0:Nx-1)+cp*  pko(0:Nx-1)
        hxk(0:Nx-1)=cx*hxk(0:Nx-1)+cp*  hpk(0:Nx-1)

        xkHxk=sum(conjg(xk(:))*hxk(:))*dx
        write(8,'(<3>i,<2>e)') ik,ib,iter,abs(hav_old-xkHxk),abs(ev-xkHxk)
        if(abs(hav_old-xkHxk)<delta_cg) exit
      enddo

        write(*,'(<3>i,<2>e)') ik,ib,iter,abs(hav_old-xkHxk),abs(ev-xkHxk)

      s=1.0d0/sqrt(sum(abs(xk(:))**2)*dx)
      u(0:Nx-1,ik,ib)=xk(0:Nx-1)*s

      u_in(0:Nx-1)=u(0:Nx-1,ik,ib)
      call h0_operation()

      xkHxk=sum(conjg(u(:,ik,ib))*hu(:))*dx
      write(*,*)
    end do
  end do
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

  do ik=-LNk,RNk
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
