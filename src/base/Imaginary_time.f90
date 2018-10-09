SUBROUTINE Imaginary_time_method()
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  integer, parameter :: Nexp_base=4
  real(8),dimension(:) :: fst(0:Nx-1)
  real(8) :: coef, rho_old, hav_old, xkxk, xkHxk
  complex(kind(0d0)),dimension(:) :: xk(0:Nx-1), hxk(0:Nx-1)
  integer :: iter, ib, jb, iexp
  integer, save :: ik
  hav_old=0.d0
  rho_old=0.d0
  do ik = 1, Nk, 1
    k_in=k(ik)
    do ib = 1, Nb, 1
      !Initial condition of wave function
      call random_number(fst)
      u(:,ik,ib)=fst(:)/sqrt(sum(abs(fst(:))**2)*dx)
      !Gram-Schmidt process
      do jb = 1, ib-1, 1
        u(:,ik,ib)=u(:,ik,ib)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*u(:,ik,ib))*dx
      end do
      xkxk=sum(conjg(u(:,ik,ib))*u(:,ik,ib))*dx
      xk(:)=u(:,ik,ib)/sqrt(xkxk)
      do iter = 1, Niter, 1
        !Imaginary time propagation!!!!!!!!!!!!!!!!
        u_in(0:Nx-1)=xk(0:Nx-1)
        coef=1.d0
        do iexp = 1, Nexp_base, 1

          coef=coef*(-1.d0)*dt/iexp
          call h0_operation()

          xk(:)=xk(:)+coef*hu(:)
          u_in(0:Nx-1)=hu(:)

        end do

        do jb = 1, ib-1, 1
          xk(:)=xk(:)-u(:,ik,jb)*sum(conjg(u(:,ik,jb))*xk(:))*dx
        end do

        xkxk =sum(conjg(xk(:))*xk(:))*dx
        xk(:)=xk(:)/sqrt(xkxk)
        u_in(0:Nx-1)=xk(:)
        call h0_operation()
        xkHxk=sum(conjg(xk(:)*hu(:))*dx)


        if((mod(iter,1000))==0)then
          write(8,'(<3>i,<3>e)') ik,ib,iter,xkxk,(xkHxk-hav_old),(abs(xk(Nx/2))**2-rho_old)
        end if
        !write(8,'(<3>i,<4>e)') ik,ib,iter,xkxk,xkHxk,(xkHxk-hav_old),(abs(xk(Nx/2))**2-rho_old)

        hav_old=xkHxk
        rho_old=abs(conjg(xk(Nx/2))*xk(Nx/2))
      end do
      u(:,ik,ib)=xk(:)
    end do
  end do

END SUBROUTINE
